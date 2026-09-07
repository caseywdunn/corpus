#!/usr/bin/env python3
"""Audit the served missing-reference ranking against raw evidence (#155).

This is a measurement tool, not another reconciliation pass. It reads a v1.3
``biblio_authority.sqlite`` and answers three separable questions:

* How many current bibliography observations exist, and how did the resolver
  map them?
* What does the frozen ``get_missing_references`` materialization rank?
* Which ranked missing works share a substantive title/year with an in-corpus
  work, and which of those also have the exact author set required to join?

The review signal is deliberately narrow: an equal normalized title and equal
year, with at least 25 alphabetic title characters. ``author_set_match`` marks
the still narrower deterministic identity rule; title/year candidates without
that evidence are not safe to merge (#225). Candidates are reported for review
and never written back to the authority database.

Every reported work includes the observation-level match methods, producer
versions and bounded source-evidence examples. That evidence explains how it
entered the ranking without changing the frozen MCP response shape.

Usage::

    python tools/qc/reference_reconciliation.py \
        --db <output>/biblio_authority.sqlite [--limit 50] [--out report.json]
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT))

from bib.authority import (  # noqa: E402
    IDENTITY_TITLE_MIN_ALPHA,
    normalize_for_key,
)


def _alpha_len(value: str) -> int:
    return sum(character.isalpha() for character in value)


def _table_names(conn: sqlite3.Connection) -> set[str]:
    return {
        row[0] for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }


def _authors_by_work(conn: sqlite3.Connection) -> dict[str, list[dict]]:
    authors: dict[str, list[dict]] = defaultdict(list)
    for row in conn.execute(
        """SELECT work_id, surname, forename, position
           FROM work_authors ORDER BY work_id, position"""
    ):
        authors[row["work_id"]].append({
            "surname": row["surname"],
            "forename": row["forename"],
            "position": row["position"],
        })
    return authors


def _author_set(authors: list[dict]) -> frozenset[str]:
    return frozenset(
        normalized
        for author in authors
        if (normalized := normalize_for_key(author.get("surname") or ""))
    )


def _current_observation_rows(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute(
        """SELECT ro.observation_id, ro.citing_corpus_hash,
                  ro.raw_citation, ro.title AS observed_title, ro.year,
                  ro.journal, ro.doi, ro.authors_json,
                  ow.work_id, ow.match_method, ow.match_score,
                  ow.producer_version
           FROM reference_current_sets current
           JOIN reference_observation_memberships member
             ON member.corpus_hash = current.corpus_hash
            AND member.source_fingerprint = current.source_fingerprint
           JOIN reference_observations ro
             ON ro.observation_id = member.observation_id
           LEFT JOIN observation_work ow
             ON ow.observation_id = ro.observation_id
           ORDER BY ro.observation_id"""
    ).fetchall()


def build_report(
    db_path: Path,
    *,
    min_citations: int = 2,
    limit: Optional[int] = 50,
    examples: int = 3,
) -> dict[str, Any]:
    """Return an evidence-backed audit of the missing-reference ranking."""
    if min_citations < 1:
        raise ValueError("min_citations must be at least 1")
    if limit is not None and limit < 0:
        raise ValueError("limit must be non-negative or None")
    if examples < 0:
        raise ValueError("examples must be non-negative")
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        required = {
            "works", "work_authors", "citations", "reference_observations",
            "reference_observation_memberships", "reference_current_sets",
            "observation_work",
        }
        absent = sorted(required - _table_names(conn))
        if absent:
            raise ValueError(
                "authority database predates the v1.3 observation schema; "
                "re-run bib.authority first (missing: " + ", ".join(absent) + ")"
            )

        observations = _current_observation_rows(conn)
        observation_count = len(observations)
        mapped = [row for row in observations if row["work_id"] is not None]
        unmapped = [row for row in observations if row["work_id"] is None]
        method_counts = Counter(row["match_method"] for row in mapped)
        producer_counts = Counter(row["producer_version"] for row in mapped)
        observations_by_work: dict[str, list[sqlite3.Row]] = defaultdict(list)
        for row in mapped:
            observations_by_work[row["work_id"]].append(row)

        authors = _authors_by_work(conn)
        in_corpus_by_title_year: dict[tuple[str, Any], list[dict]] = defaultdict(list)
        for row in conn.execute(
            """SELECT work_id, title, year, journal, doi, corpus_hash
               FROM works WHERE in_corpus = 1 ORDER BY work_id"""
        ):
            title = row["title"] or ""
            normalized = normalize_for_key(title)
            if _alpha_len(normalized) < IDENTITY_TITLE_MIN_ALPHA:
                continue
            candidate = dict(row)
            candidate["authors"] = authors.get(row["work_id"], [])
            in_corpus_by_title_year[(normalized, row["year"])].append(candidate)

        missing_rows = conn.execute(
            """SELECT w.work_id, w.title, w.year, w.journal, w.doi,
                      w.guid_type, COUNT(*) AS cited_by_count
               FROM citations c JOIN works w ON c.cited_work_id = w.work_id
               WHERE w.in_corpus = 0
               GROUP BY c.cited_work_id
               HAVING COUNT(*) >= ?
               ORDER BY cited_by_count DESC, w.work_id""",
            (int(min_citations),),
        ).fetchall()

        details = []
        possible_duplicate_count = 0
        possible_duplicate_citations = 0
        exact_identity_count = 0
        exact_identity_citations = 0
        for rank, work in enumerate(missing_rows, start=1):
            work_id = work["work_id"]
            work_observations = observations_by_work.get(work_id, [])
            title = work["title"] or ""
            normalized = normalize_for_key(title)
            candidates = []
            if _alpha_len(normalized) >= IDENTITY_TITLE_MIN_ALPHA:
                missing_author_set = _author_set(authors.get(work_id, []))
                for candidate in in_corpus_by_title_year.get(
                    (normalized, work["year"]), []
                ):
                    annotated = dict(candidate)
                    annotated["author_set_match"] = bool(
                        missing_author_set
                        and missing_author_set == _author_set(candidate["authors"])
                    )
                    candidates.append(annotated)
            if candidates:
                possible_duplicate_count += 1
                possible_duplicate_citations += work["cited_by_count"]
            if any(candidate["author_set_match"] for candidate in candidates):
                exact_identity_count += 1
                exact_identity_citations += work["cited_by_count"]

            decision_counts = Counter(
                (
                    row["match_method"], row["match_score"],
                    row["producer_version"],
                )
                for row in work_observations
            )
            decisions = [
                {
                    "method": method,
                    "score": score,
                    "producer_version": producer,
                    "observation_count": count,
                }
                for (method, score, producer), count in sorted(
                    decision_counts.items(), key=lambda item: (
                        str(item[0][0]), str(item[0][2]), item[0][1] or 0,
                    )
                )
            ]
            evidence_examples = []
            seen_examples = set()
            if examples > 0:
                for row in work_observations:
                    evidence = {
                        "observation_id": row["observation_id"],
                        "raw_citation": row["raw_citation"] or "",
                        "title": row["observed_title"] or "",
                        "year": row["year"],
                        "journal": row["journal"] or "",
                        "doi": row["doi"] or "",
                        "authors": json.loads(row["authors_json"] or "[]"),
                    }
                    signature = json.dumps(
                        {key: value for key, value in evidence.items()
                         if key != "observation_id"},
                        sort_keys=True,
                        default=str,
                    )
                    if signature not in seen_examples:
                        seen_examples.add(signature)
                        evidence_examples.append(evidence)
                    if len(evidence_examples) >= examples:
                        break
            details.append({
                "rank": rank,
                **dict(work),
                "authors": authors.get(work_id, []),
                "observation_count": len(work_observations),
                "citing_paper_count": len({
                    row["citing_corpus_hash"] for row in work_observations
                }),
                "mapping_decisions": decisions,
                "observation_examples": evidence_examples,
                "possible_in_corpus_matches": candidates,
                "review_reason": (
                    (
                        "same_normalized_title_year_and_author_set"
                        if any(c["author_set_match"] for c in candidates)
                        else "same_normalized_title_and_year"
                    )
                    if candidates else None
                ),
            })

        shown = details if limit is None else details[:max(0, int(limit))]
        compatibility_edges = conn.execute(
            "SELECT COUNT(*) FROM citations"
        ).fetchone()[0]
        return {
            "schema_version": 1,
            "database": str(db_path),
            "parameters": {
                "min_citations": int(min_citations),
                "limit": limit,
                "examples": int(examples),
                "possible_duplicate_rule": (
                    "same normalized title and year; title has at least "
                    f"{IDENTITY_TITLE_MIN_ALPHA} alphabetic characters; "
                    "author_set_match marks the resolver-safe subset"
                ),
            },
            "population": {
                "current_reference_observations": observation_count,
                "mapped_current_observations": len(mapped),
                "unmapped_current_observations": len(unmapped),
                "unmapped_citing_documents": len({
                    row["citing_corpus_hash"] for row in unmapped
                }),
                "compatibility_citation_edges": compatibility_edges,
                "observations_collapsed_by_compatibility_key": (
                    len(mapped) - compatibility_edges
                ),
                "missing_works_at_threshold": len(missing_rows),
                "missing_work_citation_edges": sum(
                    row["cited_by_count"] for row in missing_rows
                ),
                "possible_in_corpus_duplicates": possible_duplicate_count,
                "possible_duplicate_citation_edges": possible_duplicate_citations,
                "exact_identity_candidates": exact_identity_count,
                "exact_identity_citation_edges": exact_identity_citations,
            },
            "mapping_methods": dict(sorted(method_counts.items())),
            "producer_versions": dict(sorted(producer_counts.items())),
            "top_missing": shown,
        }
    finally:
        conn.close()


def _print_summary(report: dict[str, Any]) -> None:
    population = report["population"]
    print(
        "Current observations: "
        f"{population['mapped_current_observations']:,} mapped / "
        f"{population['current_reference_observations']:,} total; "
        f"{population['unmapped_citing_documents']:,} citing documents unmapped"
    )
    print(
        "Frozen citation edges: "
        f"{population['compatibility_citation_edges']:,} "
        "(" + str(population["observations_collapsed_by_compatibility_key"]) +
        " repeated same-paper observations materialized once)"
    )
    print(
        f"Missing works (threshold): {population['missing_works_at_threshold']:,}; "
        "same-title/year possible in-corpus duplicates: "
        f"{population['possible_in_corpus_duplicates']:,}; "
        "exact author-set identities still unresolved: "
        f"{population['exact_identity_candidates']:,}"
    )
    print("Mapping methods: " + json.dumps(report["mapping_methods"], sort_keys=True))
    for item in report["top_missing"]:
        flag = " REVIEW" if item["possible_in_corpus_matches"] else ""
        print(
            f"{item['rank']:>3}. {item['cited_by_count']:>4} citing papers "
            f"{item['work_id']} — {item['title'] or '<no title>'}{flag}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--min-citations", type=int, default=2)
    parser.add_argument(
        "--limit", type=int, default=50,
        help="number of ranked works in detail; use 0 for summary only",
    )
    parser.add_argument(
        "--examples", type=int, default=3,
        help="maximum source-observation examples per ranked work",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    try:
        report = build_report(
            args.db,
            min_citations=args.min_citations,
            limit=args.limit,
            examples=args.examples,
        )
    except (OSError, sqlite3.Error, ValueError) as error:
        parser.error(str(error))

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    if not args.quiet:
        _print_summary(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

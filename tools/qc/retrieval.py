#!/usr/bin/env python3
"""Capture and source-score frozen retrieval workflows; never tune ranking (#320)."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import sys
import unicodedata


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def normalize(text):
    return " ".join(unicodedata.normalize("NFKC", text).casefold().split())


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def validate_manifest(manifest):
    queries = manifest["queries"]
    if not 1 <= len(queries) <= 200:
        raise ValueError("manifest must contain 1–200 queries")
    ids = [q["id"] for q in queries]
    if len(ids) != len(set(ids)):
        raise ValueError("query IDs must be unique")
    for q in queries:
        call = q["call"]
        if (not isinstance(call.get("query"), str) or not 1 <= call.get("k", 5) <= 100
                or call.get("with_text") is not True):
            raise ValueError("each query needs text, with_text=true and k in 1–100")
        for target in q.get("targets", []):
            if target["grade"] not in {0, 1, 2}:
                raise ValueError("source relevance grades are 0, 1 or 2")
            if target["review_status"] == "source_verified":
                if not target.get("source") or not target.get("all_text") or not target.get("paper_hash"):
                    raise ValueError("source-verified targets require source evidence and text anchors")


def expected_calls(query):
    """Keep the original audited call, plus separately measured 5/10 depths."""
    original = dict(query["call"])
    return {"original": original, "at5": dict(original, k=5), "at10": dict(original, k=10)}


def required_papers(manifest):
    """Include all labelled source papers and original paper-scoped requests."""
    papers = {t["paper_hash"] for q in manifest["queries"] for t in q.get("targets", [])
              if t.get("paper_hash")}
    papers.update(q["call"]["paper_hash"] for q in manifest["queries"] if q["call"].get("paper_hash"))
    return sorted(papers)


def _check_paper_inventory(identity, required):
    """Require inventory evidence; a digest alone cannot prove target presence."""
    identity = identity if isinstance(identity, dict) else {}
    reasons = []
    population_digest = identity.get("paper_hashes_sha256")
    count = identity.get("paper_count")
    inventory = identity.get("paper_hashes")
    if (not isinstance(population_digest, str) or len(population_digest) != 64
            or any(c not in "0123456789abcdef" for c in population_digest)):
        reasons.append("missing_or_invalid_population_digest")
    if type(count) is not int or count <= 0:
        reasons.append("missing_or_invalid_population_count")
    inventory_known = (isinstance(inventory, list)
                       and all(isinstance(p, str) and p.strip() for p in inventory))
    if not inventory_known:
        reasons.append("missing_or_invalid_paper_inventory")
    else:
        if len(set(inventory)) != len(inventory):
            reasons.append("duplicate_paper_inventory_entries")
        if len(inventory) != count:
            reasons.append("paper_inventory_count_mismatch")
        if digest(sorted(inventory)) != population_digest:
            reasons.append("paper_inventory_digest_mismatch")
    identity_complete = not reasons
    required_known = (isinstance(required, list)
                      and all(isinstance(p, str) and p.strip() for p in required))
    missing = None
    if not required_known:
        reasons.append("missing_required_paper_evidence")
    elif inventory_known:
        missing = sorted(set(required)-set(inventory))
        if missing:
            reasons.append("required_papers_absent")
    return {"status": "blocked" if reasons else "pass", "reasons": reasons,
            "identity_complete": identity_complete,
            "required_paper_hashes": sorted(set(required)) if required_known else None,
            "missing_required_paper_hashes": missing}


def check_population(identity, required):
    """Check artifact membership and the papers actually competing in search."""
    identity = identity if isinstance(identity, dict) else {}
    result = _check_paper_inventory(identity, required)
    indexed = identity.get("indexed_population")
    indexed = indexed if isinstance(indexed, dict) else {}
    index_check = _check_paper_inventory(indexed, required)
    problems = []
    if indexed.get("status") != "recorded":
        problems.append("indexed_population_unavailable")
    if (indexed.get("method") != "pipeline.embedding_state.row_census"
            or indexed.get("table_name") != "document_chunks"):
        problems.append("indexed_census_method_unknown")
    if indexed.get("version_scope") != "capture":
        problems.append("indexed_population_not_bracketed_at_capture")
    before, after = indexed.get("table_version_before"), indexed.get("table_version_after")
    if type(before) is not int or type(after) is not int or before <= 0 or after <= 0:
        problems.append("indexed_table_version_unknown")
    elif before != after:
        problems.append("indexed_table_changed_during_capture")
    counts = indexed.get("paper_row_counts")
    rows = indexed.get("row_count")
    papers = indexed.get("paper_hashes")
    if (not isinstance(counts, dict) or not isinstance(papers, list)
            or not all(isinstance(p, str) for p in papers)
            or set(counts) != set(papers)
            or any(type(n) is not int or n <= 0 for n in counts.values())
            or type(rows) is not int or rows <= 0 or rows != sum(counts.values())):
        problems.append("indexed_row_counts_invalid")
    artifact_only, orphans = None, None
    if result["identity_complete"] and index_check["identity_complete"]:
        artifact_only = sorted(set(identity["paper_hashes"])-set(papers))
        orphans = sorted(set(papers)-set(identity["paper_hashes"]))
        if orphans:
            problems.append("indexed_papers_outside_artifact_inventory")
    index_check["identity_complete"] &= not problems
    index_check["reasons"].extend(problems)
    index_check["status"] = "blocked" if index_check["reasons"] else "pass"
    index_check.update(row_count=rows, artifact_only_paper_hashes=artifact_only,
                       orphan_indexed_paper_hashes=orphans)
    result["indexed_population"] = index_check
    if index_check["status"] != "pass":
        result["status"] = "blocked"
        result["reasons"].append("indexed_population_incomplete")
    return result


def _open_index_table(output_dir):
    import lancedb
    location = Path(output_dir) / "vector_db" / "lancedb"
    if not location.is_dir():
        raise FileNotFoundError("vector index directory is absent")
    return lancedb.connect(str(location)).open_table("document_chunks")


def _finish_index_population(output_dir, population):
    """Read the latest table version without loading any vectors or model."""
    try:
        population["table_version_after"] = _open_index_table(output_dir).version
        if (population.get("status") == "recorded"
                and population["table_version_after"] != population.get("table_version_before")):
            population["status"] = "changed"
    except Exception as exc:
        population["status"] = "unavailable"
        population.setdefault("error", str(exc))
        population.setdefault("exception_type", type(exc).__name__)
    return population


def capture_index_population(output_dir):
    """One projected hash/generation census; a standalone census is diagnostic."""
    from pipeline.embedding_state import row_census
    result = {"status": "unavailable", "method": "pipeline.embedding_state.row_census",
              "table_name": "document_chunks", "version_scope": "supplementary_census"}
    try:
        table = _open_index_table(output_dir)
        result["table_version_before"] = table.version
        census = row_census(table)
        counts = {paper: sum(generations.values()) for paper, generations in census.items()}
        result["row_count"] = sum(counts.values())
        if any(not isinstance(p, str) or not p.strip() for p in counts):
            raise ValueError("indexed rows lack a nonempty paper hash")
        papers = sorted(counts)
        result.update(status="recorded", paper_hashes=papers, paper_hashes_sha256=digest(papers),
                      paper_count=len(papers), paper_row_counts=counts)
    except Exception as exc:
        result.update(error=str(exc), exception_type=type(exc).__name__)
    return _finish_index_population(output_dir, result)


def target_match(target, row):
    if target.get("review_status") != "source_verified":
        return False
    if row.get("paper_hash") != target["paper_hash"]:
        return False
    text = normalize(row.get("text", ""))
    if not all(normalize(anchor) in text for anchor in target["all_text"]):
        return False
    # Old bundles lack page provenance. Exact source-graded passage anchors
    # remain usable; known conflicting page provenance is never ignored.
    pages = row.get("source_pages", [])
    return not pages or not target.get("pages") or bool(set(pages) & set(target["pages"]))


def score_rows(query, rows, k):
    rows = rows[:k]
    targets = query.get("targets", [])
    verified = [t for t in targets if t.get("review_status") == "source_verified"]
    positives = [t for t in verified if t["grade"] >= 2]
    matches = [[t for t in verified if target_match(t, r)] for r in rows]
    grades = [max((t["grade"] for t in m), default=None) for m in matches]
    docs = Counter(r.get("paper_hash") for r in rows if r.get("paper_hash"))
    seen_tables = set()
    seen_table_text = set()
    repeated = 0
    duplicate_table_text = 0
    table_rows = 0
    for row in rows:
        keys = {(row.get("paper_hash"), ref) for ref in row.get("table_refs", [])}
        table_rows += bool(keys)
        repeated += bool(keys & seen_tables)
        seen_tables.update(keys)
        if keys:
            text_key = (row.get("paper_hash"), normalize(row.get("text", "")))
            duplicate_table_text += text_key in seen_table_text
            seen_table_text.add(text_key)
    table_known = sum(bool(r.get("table_metadata_available")) for r in rows)
    covered = {t["id"] for m in matches for t in m if t["grade"] >= 2}
    return {
        "returned": len(rows), "hit": any(g is not None and g >= 2 for g in grades) if positives else None,
        "grades": grades, "unjudged_rows": grades.count(None),
        "known_positive_targets": len(positives), "positive_targets_retrieved": len(covered),
        "known_target_coverage": len(covered) / len(positives) if positives else None,
        "unique_documents": len(docs), "document_diversity": len(docs) / len(rows) if rows else 0,
        "dominant_document_fraction": max(docs.values(), default=0) / len(rows) if rows else 0,
        "table_rows_identified": table_rows, "table_metadata_rows": table_known,
        "repeated_table_rows": repeated,
        "repeated_table_rate": repeated / len(rows) if rows and table_known == len(rows) else None,
        "repeated_table_rate_known_lower_bound": repeated / len(rows) if rows else 0,
        "duplicate_table_text_rows": duplicate_table_text,
        "duplicate_table_text_rate": duplicate_table_text / len(rows) if rows and table_known == len(rows) else None,
        "matched_by_source_anchor_without_page": sum(bool(m) and not r.get("source_pages")
                                                        for m, r in zip(matches, rows)),
    }


def evaluate(manifest, capture):
    validate_manifest(manifest)
    if capture.get("manifest_sha256") != digest(manifest):
        raise ValueError("capture was made with a different manifest; freeze/review labels before capture")
    population = check_population(capture.get("identity"), required_papers(manifest))
    queries = {q["id"]: q for q in manifest["queries"]}
    reports, grouped = [], defaultdict(list)
    seen = set()
    for run in capture["runs"]:
        query = queries[run["query_id"]]
        mode = run.get("mode", "unassisted")
        if mode not in {"unassisted", "assisted"}:
            raise ValueError("run mode must explicitly distinguish assisted recovery")
        variant = run["variant"]
        key = (query["id"], mode, variant)
        if key in seen:
            raise ValueError("duplicate run variant would change the denominator")
        seen.add(key)
        if mode == "unassisted" and run["call"] != expected_calls(query).get(variant):
            raise ValueError("unassisted query or arguments differ from the frozen manifest")
        failed = any("error" in r for r in run["rows"])
        metric = score_rows(query, [] if failed else run["rows"], run["call"].get("k", 10))
        metric["execution_error"] = failed
        report = {"query_id": query["id"], "group": query["group"], "mode": mode,
                  "variant": variant, **metric}
        reports.append(report)
        if mode == "unassisted" and variant in {"at5", "at10"}:
            grouped[(query["group"], variant)].append(report)
    gates = []
    for group, budget in manifest["acceptance"]["groups"].items():
        expected = [q for q in queries.values() if q["group"] == group]
        for variant in ("at5", "at10"):
            rows = grouped[(group, variant)]
            threshold = budget["hit" + variant[2:]]
            judged = [r for r in rows if r["hit"] is not None]
            complete = (len(rows) == len(expected) and len(judged) == len(expected)
                        and len(expected) >= budget.get("minimum_queries", 1)
                        and not any(r["execution_error"] for r in rows))
            target_papers = {t["paper_hash"] for q in expected for t in q.get("targets", [])
                             if t.get("review_status") == "source_verified" and t["grade"] >= 2}
            complete &= len(target_papers) >= budget.get("minimum_papers", 1)
            hits = sum(r["hit"] is True for r in judged)
            rate = hits / len(expected) if expected else None
            gates.append({"group": group, "variant": variant, "queries": len(expected),
                          "judged_queries": len(judged), "hits": hits, "hit_rate": rate,
                          "required_rate": threshold,
                          "status": "blocked" if not complete else "pass" if rate >= threshold else "fail"})
    missing = [{"query_id": qid, "variant": variant} for qid in queries
               for variant in ("original", "at5", "at10")
               if (qid, "unassisted", variant) not in seen]
    errors = [{"query_id": r["query_id"], "variant": r["variant"]} for r in reports
              if r["mode"] == "unassisted" and r["execution_error"]]
    return {"manifest_sha256": digest(manifest), "run_identity": capture.get("identity", {}),
            "population_check": population,
            "missing_unassisted_runs": missing,
            "failed_unassisted_runs": errors,
            "status": "blocked" if (population["status"] == "blocked" or missing or errors
                                    or any(g["status"] == "blocked" for g in gates))
            else "fail" if any(g["status"] == "fail" for g in gates) else "pass",
            "gates": gates, "queries": reports,
            "interpretation": "Selected workflow benchmark; hit rates are not a deployment-wide failure rate. "
                              "Unmatched results are unjudged, not proven irrelevant. Assisted runs never count toward gates."}


def compare_reports(reference, candidate):
    if reference["manifest_sha256"] != candidate["manifest_sha256"]:
        raise ValueError("reference and candidate must use exactly the same frozen manifest")
    before = {(g["group"], g["variant"]): g for g in reference["gates"]}
    changes = []
    for gate in candidate["gates"]:
        old = before[(gate["group"], gate["variant"])]
        delta = None if None in (old["hit_rate"], gate["hit_rate"]) else gate["hit_rate"] - old["hit_rate"]
        changes.append({"group": gate["group"], "variant": gate["variant"],
                        "reference_hit_rate": old["hit_rate"], "candidate_hit_rate": gate["hit_rate"],
                        "delta": delta})
    independent = [c for c in changes if c["group"] == "independent"]
    controls = [c for c in changes if c["group"] in {"prose_control", "historical_control"}]
    checks = {}
    for name, report in (("reference", reference), ("candidate", candidate)):
        evidence = report.get("population_check")
        required = evidence.get("required_paper_hashes") if isinstance(evidence, dict) else None
        checks[name] = check_population(report.get("run_identity"), required)
    blockers = [f"{name}_population_incomplete" for name, check in checks.items() if check["status"] != "pass"]
    if all(check["identity_complete"] for check in checks.values()):
        if any(reference["run_identity"][key] != candidate["run_identity"][key]
               for key in ("paper_hashes_sha256", "paper_count")):
            blockers.append("different_paper_populations")
    if all(check["indexed_population"]["identity_complete"] for check in checks.values()):
        if any(reference["run_identity"]["indexed_population"][key]
               != candidate["run_identity"]["indexed_population"][key]
               for key in ("paper_hashes_sha256", "paper_count")):
            blockers.append("different_indexed_paper_populations")
    if checks["reference"]["required_paper_hashes"] != checks["candidate"]["required_paper_hashes"]:
        blockers.append("different_required_paper_evidence")
    complete = (not blockers and reference["status"] in {"pass", "fail"}
                and candidate["status"] in {"pass", "fail"})
    improved = (complete and bool(independent) and all(c["delta"] is not None and c["delta"] >= 0 for c in independent + controls)
                and any(c["delta"] > 0 for c in independent) and candidate["status"] == "pass")
    return {"manifest_sha256": candidate["manifest_sha256"],
            "reference_identity": reference.get("run_identity", {}), "candidate_identity": candidate.get("run_identity", {}),
            "population_comparison": {"status": "blocked" if blockers else "pass",
                                      "reasons": blockers, **checks},
            "candidate_release_gates": candidate["status"], "independent_improvement_demonstrated": improved,
            "status": "blocked" if not complete else "pass" if improved else "fail", "changes": changes}


def enrich_rows(rows, output_dir, cache):
    for row in rows:
        paper = row.get("paper_hash")
        if not paper:
            continue
        if paper not in cache:
            path = output_dir / "documents" / paper / "chunks.json"
            artifact = load(path) if path.exists() else {}
            cache[paper] = ({c["chunk_id"]: c for c in artifact.get("chunks", [])},
                            bool(artifact.get("treatment_context_policy")))
        chunks, available = cache[paper]
        chunk = chunks.get(row.get("chunk_id"), {})
        row["source_pages"] = sorted({p.get("source_page", p.get("page"))
                                       for p in chunk.get("source_items", [])
                                       if p.get("source_page", p.get("page"))})
        row["table_refs"] = sorted({t["table_ref"] for t in chunk.get("tables", []) if t.get("table_ref")})
        row["table_metadata_available"] = available and bool(chunk) and (
            bool(chunk.get("source_items")) or "tables" in chunk)
    return rows


def capture_local(manifest, output_dir, label, role):
    validate_manifest(manifest)
    # Local read-only MCP route: same query embedder and selection as serving.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from mcpsrv import app
    from mcpsrv.indexes import CorpusIndex
    from mcpsrv.tools.chunks import get_chunks_for_topic
    index = CorpusIndex(output_dir)
    index.load()
    population = capture_index_population(output_dir)
    population["version_scope"] = "capture"
    previous = app._INDEX
    app.set_index(index)
    runs, cache = [], {}
    try:
        for query in manifest["queries"]:
            for variant, call in expected_calls(query).items():
                try:
                    rows = get_chunks_for_topic(**call)
                except Exception as exc:
                    rows = [{"error": str(exc), "exception_type": type(exc).__name__}]
                runs.append({"query_id": query["id"], "mode": "unassisted", "variant": variant,
                             "call": call, "rows": enrich_rows(rows, output_dir, cache)})
    finally:
        app.set_index(previous)
        _finish_index_population(output_dir, population)
    papers = sorted(index.papers)
    return {"manifest_sha256": digest(manifest), "identity": {
        "label": label, "role": role, "output_dir": str(output_dir.resolve()),
        "captured_at": datetime.now(timezone.utc).isoformat(),
        "bundle_manifest": index.bundle_manifest, "embedding_identity": index._embedding_identity,
        "paper_hashes_sha256": digest(papers), "paper_count": len(papers), "paper_hashes": papers,
        "indexed_population": population,
        "historical_audit_equivalence": "not_asserted; compare the recorded bundle identity explicitly",
    }, "runs": runs}


def sample_independent(manifest, output_dir, count, seed):
    """Select source-review candidates before inspecting retrieval results.

    Source units, not overlapping chunks, define the population. Audit target
    regions are excluded. Selection round-robins a seeded paper order so a
    single large treatment cannot fill the held-out set.
    """
    excluded = {(t["paper_hash"], p) for q in manifest["queries"] if q["group"] != "independent"
                for t in q.get("targets", [])
                for p in t.get("pages", [])}
    population = {}
    for path in sorted((output_dir / "documents").glob("*/chunks.json")):
        paper = path.parent.name
        for row in load(path).get("chunks", []):
            treatment = row.get("treatment_context", {})
            is_key = bool(row.get("key_branches")) or any(t.get("kind") == "identification_key" for t in row.get("tables", []))
            kind = "key" if is_key else "diagnosis" if row.get("section_type") == "diagnosis" else None
            pages = sorted({i.get("source_page", i.get("page")) for i in row.get("source_items", [])
                            if i.get("source_page", i.get("page"))})
            if not kind or not pages or any((paper, p) in excluded for p in pages):
                continue
            subject = treatment.get("name") if treatment.get("status") == "resolved" else None
            if not subject and is_key:
                subject = " / ".join(row.get("headings", []))
            if not subject:
                continue
            unit = (paper, kind, subject, tuple(pages))
            population.setdefault(unit, {"paper_hash": paper, "kind": kind, "subject": subject,
                                         "pages": pages, "example_chunk_id": row["chunk_id"],
                                         "source_items": row.get("source_items", [])})
    pool = [population[key] for key in sorted(population)]
    rng = random.Random(seed)
    by_paper = defaultdict(list)
    for candidate in pool:
        by_paper[candidate["paper_hash"]].append(candidate)
    papers = sorted(by_paper)
    rng.shuffle(papers)
    for candidates in by_paper.values():
        rng.shuffle(candidates)
    selected = []
    while papers and len(selected) < count:
        for paper in list(papers):
            if len(selected) >= count:
                break
            selected.append(by_paper[paper].pop())
            if not by_paper[paper]:
                papers.remove(paper)
    result = json.loads(json.dumps(manifest))
    result["queries"] = [q for q in result["queries"] if q["group"] != "independent"]
    for n, item in enumerate(selected, 1):
        query = f"{item['subject']} diagnosis" if item["kind"] == "diagnosis" else f"{item['subject']} identification key distinguishing characters"
        result["queries"].append({"id": f"independent_{n:02d}", "group": "independent",
                                  "call": {"query": query, "k": 5, "with_text": True},
                                  "source_review_candidate": item, "targets": []})
    result["independent_sampling"] = {"seed": seed, "population_sha256": digest(pool),
        "eligible_units": len(pool), "eligible_papers": len(by_paper), "requested": count,
        "selected": len(selected), "selected_papers": len({s["paper_hash"] for s in selected}),
        "rule": "Materialized diagnosis/key source units; exclude fixed-query target pages; seeded paper round-robin.",
        "review_status": "pending_source_review; freeze source-graded targets before any ranking experiment",
        "selection": selected}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("capture", "score", "sample", "compare"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--out", required=True, type=Path)
        if name == "compare":
            cmd.add_argument("--reference", required=True, type=Path)
            cmd.add_argument("--candidate", required=True, type=Path)
            continue
        cmd.add_argument("--manifest", required=True, type=Path)
        if name != "score":
            cmd.add_argument("--output-dir", required=True, type=Path)
        if name == "capture":
            cmd.add_argument("--label", required=True)
            cmd.add_argument("--role", choices=("reference", "candidate"), required=True)
        if name == "score":
            cmd.add_argument("--capture", required=True, type=Path)
        if name == "sample":
            cmd.add_argument("--count", type=int, default=20)
            cmd.add_argument("--seed", type=int, default=3202026)
    args = parser.parse_args()
    if args.command == "compare":
        result = compare_reports(load(args.reference), load(args.candidate))
        save(args.out, result)
        return 0 if result["status"] == "pass" else 2
    manifest = load(args.manifest)
    validate_manifest(manifest)
    if args.command == "capture":
        result = capture_local(manifest, args.output_dir, args.label, args.role)
    elif args.command == "score":
        result = evaluate(manifest, load(args.capture))
    else:
        if not 1 <= args.count <= 100:
            parser.error("--count must be in 1–100")
        result = sample_independent(manifest, args.output_dir, args.count, args.seed)
    save(args.out, result)
    if args.command == "score" and result["status"] != "pass":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

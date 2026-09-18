"""Adjudicate publication years using independent source evidence (#314).

A date embedded in a title is not a publication year. Cross-year candidates
require a complete author-set match and exact title agreement after removing
that specific title date; reassignment additionally needs the publication year
in the raw reference's author prefix and the full canonical title in its body.
"""
from __future__ import annotations

import re
from collections import defaultdict


def candidate_index(conn):
    from .authority import normalize_for_key
    authors = defaultdict(set)
    for wid, surname in conn.execute("SELECT work_id,surname_normalized FROM work_authors"):
        if surname:
            authors[wid].add(surname)
    result = defaultdict(list)
    for wid, title, year, doi in conn.execute("""SELECT work_id,title,year,doi FROM works
            WHERE in_corpus=1 AND bib_imported_at IS NOT NULL
              AND title IS NOT NULL AND year IS NOT NULL ORDER BY work_id"""):
        result[frozenset(authors[wid])].append({"work_id": wid, "title": title,
            "normalized_title": normalize_for_key(title), "year": year, "doi": doi})
    return result


def adjudicate(ref, index):
    from .authority import _normalized_ref_author_set, normalize_doi, normalize_for_key
    from .reference_quality import author_quality_reasons
    if author_quality_reasons(ref):
        return None, []
    parsed_year = ref.get("year")
    title = normalize_for_key(ref.get("title") or "")
    authors = _normalized_ref_author_set(ref.get("authors") or [])
    if not isinstance(parsed_year, int) or not authors or sum(c.isalpha() for c in title) < 25:
        return None, []
    title_date = re.compile(rf"\b{parsed_year}\b")
    stripped_title = " ".join(title_date.sub("", title).split())
    candidates = []
    for candidate in index.get(authors, []):
        if candidate["year"] == parsed_year or not title_date.search(candidate["normalized_title"]):
            continue
        stripped_candidate = " ".join(title_date.sub("", candidate["normalized_title"]).split())
        if stripped_title != stripped_candidate:
            continue
        if ref.get("doi") and normalize_doi(ref["doi"]) != normalize_doi(candidate["doi"] or ""):
            continue
        candidates.append(candidate)
    if not candidates:
        return None, []
    raw = normalize_for_key(ref.get("raw") or "")
    supported = []
    for candidate in candidates:
        start = raw.find(candidate["normalized_title"])
        if start < 0 or start > 500:
            continue
        prefix = raw[:start]
        years = {int(y) for y in re.findall(r"\b(?:1[5-9]\d{2}|20\d{2})\b", prefix)}
        if years != {candidate["year"]}:
            continue
        if not all(re.search(rf"\b{re.escape(surname)}\b", prefix) for surname in authors):
            continue
        supported.append(candidate)
    if len(supported) == 1 and "#part:" not in supported[0]["work_id"]:
        candidate = supported[0]
        return candidate["work_id"], [{"code": "raw_publication_year_supported",
            "parsed_year": parsed_year, "publication_year": candidate["year"],
            "candidate_work_id": candidate["work_id"],
            "basis": "complete_authors_and_full_canonical_title_with_publication_year_in_raw_author_prefix"}]
    return None, [{"code": "possible_publication_year_conflict", "parsed_year": parsed_year,
        "candidate_work_ids": [c["work_id"] for c in candidates[:5]],
        "candidate_count": len(candidates), "requires_source_review": True,
        "basis": "parsed_year_occurs_in_a_curated_same_author_title_but_raw_publication_evidence_is_insufficient_or_ambiguous"}]

"""Adjudicate publication years using independent source evidence (#314).

A date embedded in a title is not a publication year. Cross-year candidates
require a complete author-set match and title evidence after removing that
specific title date. Reassignment needs the publication year in the raw
reference's author prefix. A single omitted article additionally requires the
canonical volume and page range in that same raw citation.
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
    for wid, title, year, doi, volume, pages in conn.execute("""SELECT work_id,title,year,doi,volume,pages FROM works
            WHERE in_corpus=1 AND bib_imported_at IS NOT NULL
              AND title IS NOT NULL AND year IS NOT NULL ORDER BY work_id"""):
        result[frozenset(authors[wid])].append({"work_id": wid, "title": title,
            "normalized_title": normalize_for_key(title), "year": year, "doi": doi,
            "volume": volume, "pages": pages})
    return result


def _title_forms(title):
    """Exact title, then one omitted 'the'; no substantive-word fuzziness."""
    yield title, "exact"
    words = title.split()
    for i, word in enumerate(words):
        if word == "the":
            yield " ".join(words[:i] + words[i + 1:]), "single_article_omission"


def _locator_supported(raw_tail, candidate, ref):
    """Require a contiguous volume/page pair, never numbers anywhere in text."""
    volume = str(candidate.get("volume") or "").strip()
    pages = re.fullmatch(r"(\d+)\s*[-–—]+\s*(\d+)", str(candidate.get("pages") or ""))
    if not volume.isdigit() or not pages:
        return False
    first, last = pages.groups()
    # Some producers expose parsed locators as well as raw text. Contradictory
    # fields need review; this rule is not a general locator-repair mechanism.
    if ref.get("volume") and str(ref["volume"]).strip() != volume:
        return False
    if ref.get("pages"):
        parsed_pages = re.fullmatch(r"(\d+)\s*[-–—]+\s*(\d+)", str(ref["pages"]).strip())
        if parsed_pages is None or parsed_pages.groups() != (first, last):
            return False
    return bool(re.search(rf"(?<!\w){re.escape(volume)}\s*[,;:]\s*"
                          rf"{re.escape(first)}\s*[-–—]+\s*{re.escape(last)}"
                          rf"(?!\w|\s*[-–—])", raw_tail))


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
        if stripped_title not in {form for form, _ in _title_forms(stripped_candidate)}:
            continue
        if ref.get("doi") and normalize_doi(ref["doi"]) != normalize_doi(candidate["doi"] or ""):
            continue
        candidates.append(candidate)
    if not candidates:
        return None, []
    raw_source = ref.get("raw") or ""
    # Grobid flattens a printed line break to a space: ``siphono- phores``.
    # This comparison-only join still has to reproduce the entire curated
    # title; neither the immutable observation nor the source prose is edited.
    title_words = {word for candidate in candidates for word in candidate["normalized_title"].split()}
    def join_title_word(match):
        joined = match[1] + match[2]
        return joined if normalize_for_key(joined) in title_words else match[0]
    raw_joined = re.sub(r"\b([^\W\d_]+)-\s+([^\W\d_]+)\b", join_title_word, raw_source)
    raw = normalize_for_key(raw_joined)
    supported = []
    for candidate in candidates:
        for form, alignment in _title_forms(candidate["normalized_title"]):
            match = re.search(rf"(?<!\w){re.escape(form)}(?!\w)", raw)
            if match is None or match.start() > 500:
                continue
            prefix = raw[:match.start()]
            years = {int(y) for y in re.findall(r"\b(?:1[5-9]\d{2}|20\d{2})\b", prefix)}
            if years != {candidate["year"]}:
                continue
            if not all(re.search(rf"\b{re.escape(surname)}\b", prefix) for surname in authors):
                continue
            if alignment != "exact":
                # The prefix has already excluded this date, so its first
                # occurrence ends the matched title. Preserve locator dashes.
                title_date_end = re.search(rf"\b{parsed_year}\b", raw_joined)
                if title_date_end is None or not _locator_supported(raw_joined[title_date_end.end():], candidate, ref):
                    continue
            supported.append((candidate, alignment))
            break
    if len(supported) == 1 and "#part:" not in supported[0][0]["work_id"]:
        candidate, alignment = supported[0]
        return candidate["work_id"], [{"code": "raw_publication_year_supported",
            "parsed_year": parsed_year, "publication_year": candidate["year"],
            "candidate_work_id": candidate["work_id"],
            "raw_title_alignment": alignment,
            "line_hyphen_joined_for_comparison": raw_joined != raw_source,
            "basis": ("complete_authors_and_full_canonical_title_with_publication_year_in_raw_author_prefix"
                      if alignment == "exact" else
                      "complete_authors_and_single_article_omission_with_publication_year_and_volume_pages_in_raw")}]
    return None, [{"code": "possible_publication_year_conflict", "parsed_year": parsed_year,
        "candidate_work_ids": [c["work_id"] for c in candidates[:5]],
        "candidate_count": len(candidates), "requires_source_review": True,
        "basis": "parsed_year_occurs_in_a_curated_same_author_title_but_raw_publication_evidence_is_insufficient_or_ambiguous"}]

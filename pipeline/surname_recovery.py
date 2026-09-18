"""Conservative, source-local recovery of citation surnames (#315).

Curated author/year pairs only propose candidates. Unhinted regional OCR must
agree on the complete surname in two segmentation modes, with an independently
valid matching year and no conflicting date. No name is injected into OCR
dictionaries. Printed alternatives remain untouched.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import unicodedata

SURNAME_POLICY = "curated-author-year-source-ocr-consensus-v3"
_CITATION = re.compile(
    r"(?<!\w)([^\W\d_]{5,})\s*(?:,\s*|\(\s*|et\s+al\.\s*)?((?:1[5-9]|20)\d{2})[a-z]?(?!\d)",
    re.UNICODE,
)

# Read complete numeric tokens, including damaged suffixes, without treating a
# prefix of e.g. 19814 as the year 1981. The independent name and year evidence
# is adjudicated below; this pattern does not change the document's date text.
_OCR_NAME_NUMBER = re.compile(
    r"(?<!\w)([^\W\d_]{5,})\s*(?:,\s*|\(\s*|et\s+al\.\s*)?([0-9]\w*)(?!\w)",
    re.UNICODE,
)
_OCR_NUMBER = re.compile(r"(?<!\d)([0-9]\w*)(?!\w)")
_OCR_YEAR = re.compile(r"((?:1[5-9]|20)[0-9]{2})[a-z]?")


def _key(value):
    return unicodedata.normalize("NFC", value).casefold()


def author_catalog(entries):
    """Small deterministic reviewed catalog, independent of PDF membership."""
    from bib.parser import _split_authors, entry_ocrlang

    pairs = {}
    languages = defaultdict(Counter)
    language_sources = defaultdict(list)
    for entry in entries:
        year = str(entry.get("year") or "")
        if not re.fullmatch(r"(?:1[5-9]|20)\d{2}", year):
            continue
        for author in _split_authors(entry.get("author") or ""):
            surname = unicodedata.normalize("NFC", author["surname"])
            if not surname:
                continue
            key = (_key(surname), int(year))
            row = pairs.setdefault(
                key,
                {
                    "surname": surname,
                    "year": int(year),
                    "sources": [],
                    "_languages": set(),
                },
            )
            row["surname"] = min(
                (row["surname"], surname), key=lambda value: (value.isupper(), value)
            )
            row["sources"].append(
                {"bib_key": entry.get("_key"), "title": entry.get("title") or ""}
            )
            declared = set((entry_ocrlang(entry) or "eng").split("+")) - {"eng", "osd"}
            row["_languages"].update(declared)
            languages[_key(surname)].update(declared)
            if declared:
                language_sources[_key(surname)].append(
                    {
                        "bib_key": entry.get("_key"),
                        "year": int(year),
                        "languages": sorted(declared),
                    }
                )
    records = []
    for key, row in sorted(pairs.items()):
        if len(row["surname"]) < 5 or not any(ord(c) > 127 for c in row["surname"]):
            # These pairs only exclude already known names in propose(). Their
            # titles, entry keys and language declarations cannot support an
            # OCR candidate or reference repair, so do not invalidate extraction
            # for ordinary title/key curation on these records.
            records.append(
                {
                    "surname": row["surname"],
                    "year": row["year"],
                    "sources": [],
                    "languages": [],
                    "language_basis": "known_name_exclusion",
                    "language_sources": [],
                }
            )
            continue
        row["sources"] = sorted(
            row["sources"], key=lambda s: (s["bib_key"] or "", s["title"])
        )
        # A source-language pack supplied in the reviewed library is preferable
        # to English for its author's distinctive spelling. No inferred locale
        # is derived from a surname or a special-case spelling substitution.
        direct = row.pop("_languages")
        choices = direct or set(languages[key[0]]) or {"eng"}
        row["languages"] = sorted(
            choices, key=lambda lang: (-languages[key[0]][lang], lang)
        )[:1]
        row["language_basis"] = (
            "same_year_curated_publication"
            if direct
            else "same_author_curated_publications"
            if languages[key[0]]
            else "default_english"
        )
        row["language_sources"] = sorted(
            language_sources[key[0]], key=lambda s: (s["bib_key"] or "", s["year"])
        )
        records.append(row)
    canonical = json.dumps(
        records, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    )
    return {
        "records": records,
        "sha256": hashlib.sha256(canonical.encode()).hexdigest(),
    }


@lru_cache(maxsize=64)
def _hash(path, mtime, size):
    del mtime, size
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def surname_recovery_producer(catalog):
    executable = shutil.which("tesseract")
    if not executable:
        sibling = Path(sys.executable).parent / "tesseract"
        executable = str(sibling) if sibling.is_file() else None
    result = {
        "policy": SURNAME_POLICY,
        "catalog_sha256": catalog["sha256"],
        "dpi": 600,
        "ocr_modes": [6, 13],
        "dictionary_hints": False,
        "max_candidates_per_page": 12,
        "max_candidates_per_document": 64,
        "timeout_seconds": 15,
        "max_crop_pixels": 2_000_000,
        "models": {},
        "pymupdf_version": importlib.metadata.version("pymupdf"),
        "available": False,
    }
    if not executable:
        return result
    try:
        version = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        ).stdout.strip()
        listing = subprocess.run(
            [executable, "--list-langs"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        ).stdout
        location = re.search(r'"([^"]+)"', listing)
        for lang in sorted(
            {lang for row in catalog["records"] for lang in row["languages"]}
        ):
            model = Path(location[1]) / f"{lang}.traineddata" if location else None
            if model and model.is_file():
                stat = model.stat()
                result["models"][lang] = _hash(
                    str(model), stat.st_mtime_ns, stat.st_size
                )
        result.update(
            available=bool(result["models"]), executable=executable, version=version
        )
    except (OSError, subprocess.SubprocessError):
        pass
    return result


def _distance(left, right):
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        current = [i]
        for j, b in enumerate(right, 1):
            current.append(
                min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (a != b))
            )
        previous = current
    return previous[-1]


def propose(text, catalog):
    """Citation-shaped near names are review leads, never corrections alone."""
    by_year = defaultdict(list)
    for row in catalog["records"]:
        by_year[row["year"]].append(row)
    proposals = []
    for match in _CITATION.finditer(text):
        observed, year = match.group(1), int(match.group(2))
        known = by_year[year]
        if any(_key(row["surname"]) == _key(observed) for row in known):
            continue
        candidates = [
            row
            for row in known
            if len(row["surname"]) >= 5
            and _key(row["surname"])[:3] == _key(observed)[:3]
            and any(ord(c) > 127 for c in row["surname"])
            and 0 < _distance(_key(observed), _key(row["surname"])) <= 2
        ]
        if candidates:
            context = text[max(0, match.start() - 2) : min(len(text), match.end() + 15)]
            quoted = bool(re.search(r'["“”«»]|\[sic\]', context, re.IGNORECASE))
            proposals.append(
                {
                    "charspan": list(match.span(1)),
                    "citation_span": list(match.span()),
                    "original": observed,
                    "year": year,
                    "candidates": candidates,
                    "status": "quoted_or_sic_context" if quoted else "candidate",
                }
            )
    return proposals


def _source_crop(page, bbox, proposal):
    """Require a unique native word/year anchor within this exact item."""
    import fitz

    region = fitz.Rect(bbox)
    words = page.get_text("words")
    hits = []
    for i, word in enumerate(words):
        token = word[4].strip("()[]{}.,;:")
        displayed = fitz.Rect(word[:4]) * page.rotation_matrix
        if _key(token) != _key(proposal["original"]) or not region.contains(
            fitz.Point(
                (displayed.x0 + displayed.x1) / 2, (displayed.y0 + displayed.y1) / 2
            )
        ):
            continue
        for following in words[i + 1 : i + 4]:
            if following[5:7] != word[5:7]:
                break
            if re.fullmatch(str(proposal["year"]) + r"[a-z]?[),.;:\]]*", following[4]):
                rect = (fitz.Rect(word[:4]) | fitz.Rect(following[:4])) + (
                    -2,
                    -1.5,
                    2,
                    1.5,
                )
                hits.append((rect * page.rotation_matrix) & page.rect)
                break
    return hits[0] if len(hits) == 1 else None


def _render_source_crop(page, rect, dpi):
    """Render upright source letters from a displayed-page crop rectangle."""
    rotation = page.rotation
    native_rect = rect * page.derotation_matrix
    try:
        # PyMuPDF text coordinates are unrotated; Docling/page crop coordinates
        # are displayed coordinates. Change only this in-memory view and always
        # restore it, so OCR sees the original baseline at every PDF rotation.
        if rotation:
            page.set_rotation(0)
        return page.get_pixmap(clip=native_rect, dpi=dpi, alpha=False).tobytes("png")
    finally:
        if rotation:
            page.set_rotation(rotation)


def _ocr_reading(text):
    """Keep adjacent name/number evidence separate from valid date readings."""
    pairs = _OCR_NAME_NUMBER.findall(text)
    tokens = _OCR_NUMBER.findall(text)

    def valid_year(token):
        match = _OCR_YEAR.fullmatch(token)
        return int(match[1]) if match else None

    return {
        "names": [name for name, _ in pairs],
        "year_tokens": tokens,
        "valid_years": [
            year for token in tokens if (year := valid_year(token)) is not None
        ],
        "adjacent_years": [
            year for _, token in pairs if (year := valid_year(token)) is not None
        ],
        "invalid_year_tokens": [token for token in tokens if valid_year(token) is None],
    }


def adjudicate_crop(png, proposal, producer):
    """No user words/whitelist: preserve source-printed alternative spellings."""
    candidate = proposal["candidates"][0]
    outputs = []
    langs = [lang for lang in candidate["languages"] if lang in producer["models"]]
    if not langs:
        return {"status": "ocr_language_unavailable", "ocr_outputs": []}
    for lang in langs:
        for mode in producer["ocr_modes"]:
            try:
                run = subprocess.run(
                    [
                        producer["executable"],
                        "stdin",
                        "stdout",
                        "-l",
                        lang,
                        "--psm",
                        str(mode),
                        "-c",
                        "load_system_dawg=0",
                        "-c",
                        "load_freq_dawg=0",
                    ],
                    input=png,
                    capture_output=True,
                    timeout=producer["timeout_seconds"],
                    check=True,
                )
                text = run.stdout.decode("utf-8").strip()
            except (OSError, subprocess.SubprocessError) as exc:
                return {
                    "status": "ocr_failed",
                    "error": type(exc).__name__,
                    "ocr_outputs": outputs,
                }
            outputs.append(
                {"language": lang, "psm": mode, "text": text, **_ocr_reading(text)}
            )
    expected = _key(candidate["surname"])
    year_evidence = {
        "source_anchor_year": proposal["year"],
        "matching_readings": [
            {"language": row["language"], "psm": row["psm"]}
            for row in outputs
            if proposal["year"] in row["adjacent_years"]
        ],
        "conflicting_valid_years": sorted(
            {
                year
                for row in outputs
                for year in row["valid_years"]
                if year != proposal["year"]
            }
        ),
        "invalid_tokens": [
            token for row in outputs for token in row["invalid_year_tokens"]
        ],
        "basis": "unique_native_anchor_and_at_least_one_complete_matching_ocr_year",
    }
    supported_year = (
        bool(year_evidence["matching_readings"])
        and not year_evidence["conflicting_valid_years"]
    )
    distinct_modes = len({row["psm"] for row in outputs}) >= 2
    result = {"ocr_outputs": outputs, "year_evidence": year_evidence}
    if (
        distinct_modes
        and supported_year
        and all(
            len(row["names"]) == 1 and _key(row["names"][0]) == expected
            for row in outputs
        )
    ):
        return {"status": "verified", "replacement": candidate["surname"], **result}
    if (
        distinct_modes
        and supported_year
        and all(
            len(row["names"]) == 1
            and _key(row["names"][0]) == _key(proposal["original"])
            for row in outputs
        )
    ):
        return {"status": "source_supports_observed_spelling", **result}
    return {"status": "ocr_disagreement", **result}


def recover_citation_surnames(document, pdf_path, catalog, *, producer=None):
    """Return persisted source decisions; mutate only corroborated text spans."""
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta

    producer = producer or surname_recovery_producer(catalog)
    source = Path(pdf_path)
    stat = source.stat()
    report = {
        "policy": SURNAME_POLICY,
        "producer": producer,
        "decisions": [],
        "catalog_sha256": catalog["sha256"],
        "source_pdf_sha256": _hash(str(source), stat.st_mtime_ns, stat.st_size),
    }
    counts = defaultdict(int)
    with fitz.open(pdf_path) as pdf:
        for item in document.texts:
            proposals = propose(item.text, catalog)
            if not proposals:
                continue
            edits = []
            for proposal in proposals:
                decision = {
                    **proposal,
                    "item_ref": item.self_ref,
                    "policy": SURNAME_POLICY,
                }
                if decision["status"] != "candidate":
                    pass
                elif len(proposal["candidates"]) != 1:
                    decision["status"] = "ambiguous_curated_names"
                elif len(item.prov) != 1 or not 1 <= item.prov[0].page_no <= len(pdf):
                    decision["status"] = "source_region_unavailable"
                else:
                    prov = item.prov[0]
                    page = pdf[prov.page_no - 1]
                    decision["page"] = prov.page_no
                    if (
                        counts[prov.page_no] >= producer["max_candidates_per_page"]
                        or sum(counts.values())
                        >= producer["max_candidates_per_document"]
                    ):
                        decision["status"] = "candidate_budget_exceeded"
                    else:
                        counts[prov.page_no] += 1
                        b = prov.bbox.to_top_left_origin(page.rect.height)
                        rect = _source_crop(page, (b.l, b.t, b.r, b.b), proposal)
                        if rect is None:
                            decision["status"] = "source_anchor_ambiguous_or_missing"
                        elif (
                            rect.width * rect.height * (producer["dpi"] / 72) ** 2
                            > producer["max_crop_pixels"]
                        ):
                            decision["status"] = "crop_pixel_budget_exceeded"
                        else:
                            png = _render_source_crop(page, rect, producer["dpi"])
                            decision.update(
                                source_bbox=list(rect),
                                source_rotation=page.rotation,
                                crop_sha256=hashlib.sha256(png).hexdigest(),
                            )
                            decision.update(adjudicate_crop(png, proposal, producer))
                if decision["status"] == "verified":
                    edits.append(decision)
                report["decisions"].append(decision)
            for edit in sorted(edits, key=lambda e: e["charspan"][0], reverse=True):
                start, end = edit["charspan"]
                replacement = edit["replacement"]
                if edit["original"].isupper():
                    replacement = replacement.upper()
                item.text = item.text[:start] + replacement + item.text[end:]
            if item.meta is None:
                item.meta = BaseMeta()
            notes = list(getattr(item.meta, "corpus__surname_recovery", []))
            for decision in report["decisions"]:
                if decision["item_ref"] == item.self_ref and decision not in notes:
                    notes.append(decision)
            item.meta.corpus__surname_recovery = notes
    report["unresolved"] = [
        decision
        for decision in report["decisions"]
        if decision["status"]
        not in {
            "verified",
            "source_supports_observed_spelling",
            "quoted_or_sic_context",
        }
    ]
    return report

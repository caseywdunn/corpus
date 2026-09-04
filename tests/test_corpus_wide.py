"""Corpus-wide content consistency checks.

These tests iterate over ALL processed documents and verify content quality
and cross-consistency — things that can be checked programmatically without
ground truth. They complement the per-paper ground-truth tests.

Categories:
  structural  — expected files exist, JSON parses, required fields present
  figures     — figure ↔ text cross-referencing consistency
  citations   — reference list ↔ corpus cross-matching
  metadata    — year/author/title plausibility via internal cross-checks
  text        — extraction quality signals (chars/page, alphabet ratio)
  chunks      — chunking quality (duplicates, empty chunks)

Run:
    python -m pytest tests/test_corpus_wide.py -v
    python -m pytest tests/test_corpus_wide.py -v -k "TestFigureTextConsistency"
    python -m pytest tests/test_corpus_wide.py -v --tb=line   # compact failures

Slow citation-graph tests (build a corpus-wide index first):
    python -m pytest tests/test_corpus_wide.py -v -k "TestCitationGraph"
"""

import json
import os
import re
from collections import defaultdict
from pathlib import Path

import pytest

pytestmark = pytest.mark.corpus_required


# ---------------------------------------------------------------------------
# Comparing metadata against extracted text
# ---------------------------------------------------------------------------


def _norm(s: str) -> str:
    """Collapse whitespace and drop space-before-punctuation for comparison.

    Metadata strings come from BibTeX or Grobid; body text comes from
    docling over a PDF (increasingly, over OCR). The two agree on words
    and disagree on spacing constantly, so comparing them literally
    tests typography rather than correctness:

    * ``Marrus claudanielis, a new species`` (bib) vs
      ``MARRUS CLAUDANIELIS , A NEW SPECIES`` (extracted) — one space
      before a comma.
    * ``Einige histologische Befunde an Coelenteraten`` vs
      ``Einige  histologische  Befunde  an  Coelenteraten.`` — doubled
      spaces from OCR.

    Both titles are plainly present. Normalizing here keeps these tests
    pointed at "is the metadata plausible for this text", which is what
    they are for, and stops them failing on any docling upgrade or
    re-OCR that shifts tokenization (#167).
    """
    s = re.sub(r"\s+", " ", (s or "").lower())
    s = re.sub(r"\s+([,;:.\)\]])", r"\1", s)
    return s.strip()


# Unicode ranges that tell us which writing system a string is in. Only
# what this corpus actually contains.
_SCRIPT_RANGES = (
    ("cyrillic", 0x0400, 0x052F),
    ("greek", 0x0370, 0x03FF),
    ("cjk", 0x3040, 0x9FFF),
)


def _script_shares(s: str) -> dict:
    """Fraction of alphabetic characters in each writing system."""
    counts = {"latin": 0}
    for ch in s or "":
        if not ch.isalpha():
            continue
        cp = ord(ch)
        for name, lo, hi in _SCRIPT_RANGES:
            if lo <= cp <= hi:
                counts[name] = counts.get(name, 0) + 1
                break
        else:
            if cp < 0x0250 or 0x1E00 <= cp <= 0x1EFF:
                counts["latin"] += 1
    total = sum(counts.values())
    if not total:
        return {}
    return {k: v / total for k, v in counts.items() if v}


def _dominant_script(s: str) -> str:
    """``'latin'``, ``'cyrillic'``, ``'greek'``, ``'cjk'`` or ``'none'``."""
    shares = _script_shares(s)
    return max(shares, key=shares.get) if shares else "none"


# A body with at least this share of non-Latin characters is a
# foreign-script paper for our purposes, whatever the raw majority says.
_NONLATIN_SUBSTANTIAL = 0.20


def _nonlatin_summary(s: str) -> str:
    """e.g. ``'42% cyrillic'`` — the non-Latin content, for skip messages."""
    shares = {k: v for k, v in _script_shares(s).items() if k != "latin"}
    if not shares:
        return "0% non-Latin"
    return ", ".join(
        f"{v:.0%} {k}" for k, v in sorted(shares.items(), key=lambda kv: -kv[1])
    )


def _substantially_nonlatin(s: str) -> bool:
    """True when enough of ``s`` is non-Latin to call it foreign-script."""
    shares = _script_shares(s)
    return sum(v for k, v in shares.items() if k != "latin") >= _NONLATIN_SUBSTANTIAL


def _transliteration_likely(meta_value: str, body: str) -> bool:
    """True when Latin metadata is describing non-Latin-script content.

    Bibliographic metadata for foreign-language papers is conventionally
    Latin — a translated title, a transliterated surname — while the body
    stays in the original. Stepanjants 1970 is recorded as "Siphonophora
    of the southern part of the Kurile-Kamchatka Trench…" and its text
    begins "СИФОНОФОРЫ РАЙОНА…"; the author is ``Stepanjants`` in the bib
    and ``Степаньянц`` on the page. Neither is wrong, and no amount of
    extraction work will make one contain the other.

    Deliberately a *share* test rather than "which script wins". Russian
    taxonomic papers are dense with Latin binomials and journal titles —
    Stepanjants 1970 is 21,084 Latin characters against 15,075 Cyrillic,
    so a majority test calls it a Latin document and the Cyrillic body
    text goes unnoticed. 42% Cyrillic is plainly a Russian paper.
    """
    if _dominant_script(meta_value) != "latin":
        return False
    return _substantially_nonlatin(body)


# ---------------------------------------------------------------------------
# Discovery and loading
# ---------------------------------------------------------------------------

# Default output location for the corpus-wide tests. Precedence
# matches tests/conftest.py — see #59.
_BOUCHET_FALLBACK = Path("/nfs/roberts/project/pi_cwd7/cwd7/output")


def _output_dir():
    override = os.environ.get("CORPUS_OUTPUT_DIR")
    if override:
        return Path(override)
    demo_output = Path(__file__).parent.parent / "demo" / "output"
    if demo_output.is_dir():
        return demo_output
    return _BOUCHET_FALLBACK


def _discover_documents():
    """Return sorted list of (hash, doc_dir) for every processed document."""
    docs_root = _output_dir() / "documents"
    if not docs_root.is_dir():
        return []
    results = []
    for d in sorted(docs_root.iterdir()):
        if d.is_dir() and (d / "summary.json").exists():
            results.append((d.name, d))
    return results


_DOCS = _discover_documents()
_DOC_IDS = [h for h, _ in _DOCS]
_DOC_MAP = dict(_DOCS)


def _doc_dir(request):
    return _DOC_MAP[request.param]


def _load_json_safe(path):
    """Load JSON, returning (data, error). Never raises."""
    try:
        with open(path) as f:
            return json.load(f), None
    except (json.JSONDecodeError, OSError) as e:
        return None, str(e)


def _load_json_or_skip(doc_dir, filename):
    path = doc_dir / filename
    if not path.exists():
        pytest.skip(f"{filename} missing")
    data, err = _load_json_safe(path)
    if err:
        pytest.fail(f"{filename} parse error: {err}")
    return data


# Pattern for figure/plate references in running text
_FIG_REF_PATTERN = re.compile(
    r"""
    (?:Fig(?:ure|s?\.?)|Plate)\s*       # "Fig.", "Figure", "Figs.", "Plate"
    (\d+)                                # first number
    (?:\s*[-–,&]\s*(\d+))*              # optional additional numbers
    """,
    re.IGNORECASE | re.VERBOSE,
)


def _extract_text_figure_numbers(text):
    """Return set of digit strings for figure/plate numbers mentioned in text.

    Returns strings (not ints) so they can be compared directly against the
    string figure_number values stored in figures.json.
    """
    nums = set()
    for m in _FIG_REF_PATTERN.finditer(text):
        # Capture range endpoints (e.g. "Figs. 1-3") as digit strings
        for n in re.findall(r"\d+", m.group(0)):
            nums.add(n)
    return nums


def _fig_number_digit_keys(fn):
    """Extract digit strings from a figure_number for comparison against text.

    Handles plain integers ('5'), subfigure panels ('3c'), and
    book-style chapter.figure numbers ('2.1').  Returns the set of
    digit-run strings extracted, e.g. '3c' -> {'3'}, '2.1' -> {'2', '1'}.
    """
    return set(re.findall(r"\d+", str(fn)))


# ---------------------------------------------------------------------------
# Corpus-wide index for citation graph tests (built once, lazily)
# ---------------------------------------------------------------------------

_CORPUS_INDEX = None

# Below this many documents, "cites nothing else in the corpus" carries no
# signal about reference parsing — it is the expected outcome. The bundled
# demo is 4 papers and none of them cites another (Pugh 2001 and Dunn 2005
# both cite *Pugh 1975/1989*, while the corpus holds *Pugh 2001*), so the
# check fired on every clean local build. A production corpuscle is
# ~1,800 papers; anything between the two is an arbitrary line, so this is
# set well above a toy fixture and far below any corpus the check is
# meant to police.
_CITATION_GRAPH_MIN_DOCS = 25


def _build_corpus_index():
    """Build a lookup of (normalized_first_author_surname, year) -> [doc_hash].

    Used to check whether references in one paper match other papers in the
    corpus. This is a fuzzy match — it won't catch everything, but
    disagreements are informative.
    """
    global _CORPUS_INDEX
    if _CORPUS_INDEX is not None:
        return _CORPUS_INDEX
    _CORPUS_INDEX = defaultdict(list)
    for doc_hash, doc_dir in _DOCS:
        md_path = doc_dir / "metadata.json"
        if not md_path.exists():
            continue
        data, _ = _load_json_safe(md_path)
        if data is None:
            continue
        year = data.get("year")
        authors = data.get("authors", [])
        if not authors or year is None:
            continue
        first = authors[0]
        surname = first.get("surname", "") if isinstance(first, dict) else str(first)
        surname = surname.strip().lower()
        if surname:
            _CORPUS_INDEX[(surname, year)].append(doc_hash)
    return _CORPUS_INDEX


# ---------------------------------------------------------------------------
# Soft consistency checks as corpus rates (#185)
# ---------------------------------------------------------------------------
#
# These checks compare two derived artifacts and read a disagreement as a
# pipeline defect. That premise holds for some of this material and not for
# the rest, and which one you get is a property of the *corpus*, not of the
# code: a curated .bib title legitimately differs from what an offprint's
# missing title page prints, and a 19th-century monograph legitimately
# references plates bound in another volume.
#
# Asserted per paper on one full production build, they produced failures
# across 65% of the corpus. A signal that
# fires on two thirds of the corpus cannot be triaged, so in practice it was
# off. Asserting the *rate* keeps the signal: 3% of born-digital papers
# failing the title check is actionable and a jump to 30% is a real
# regression, while `a1b2c3 failed` is neither.
#
# This is the same premise-drift 9a6d4b0 fixed at the small-corpus end,
# arriving from the other direction. There it was fixed by gating on corpus
# size; here corpus size is not what makes the premise wrong, so the
# assertion changes shape instead.
#
# The hard per-paper checks are deliberately left alone — has_text,
# has_title, has_authors, has_chunks, text_min_length, chars_per_page. They
# agree with the quality gates run.log already reports (low_text_density 53,
# empty_text 27) and they are few enough to act on individually.

_SOFT_RATE_MIN_DOCS = _CITATION_GRAPH_MIN_DOCS

# Ceilings are empirical, not principled. They answer "has this got
# materially worse", not "is this good", and they are calibrated against
# the two corpora available to measure — deliberately unalike, because a
# ceiling fitted to one corpus is indistinguishable from a ceiling fitted
# to that corpus's *subject matter*:
#
#   siphonophore  ~2,000 docs, marine invertebrate zoology, heavy on
#                 19th-century plate-based monographs and offprints
#   viburnum        699 docs, botany, mostly modern journal articles,
#                   assembled by mining references
#
# Observed rate (born-digital / OCR), and the ceiling set from it:
#
#   check                             siphonophore    viburnum      ceiling
#   title_appears_in_text             15.3% / 25.5%   9.3% / 19.8%   30 / 45
#   text_figure_refs_have_figures     34.0% / 66.2%   6.6% / 18.8%   50 / 80
#   first_author_surname_in_text       5.4% /  4.8%  11.6% /  7.5%   25 / 25
#   references_match_corpus_papers    26.7% / 47.8%  52.7% / 84.4%   90 / 95
#   alphabet_ratio                     4.1% /  9.0%   3.5% / 10.7%   12 / 20
#   extracted_figures_referenced       1.7% /  8.2%  11.7% /  9.1%   25 / 25
#   no_duplicate_chunks                1.0% /  0.6%   0.7% /  0.0%    5 /  5
#
# Ceiling ~= 1.5-2x the higher of the two observed rates, rounded: normal
# variation between corpora of different character must not fire, while a
# doubling of a rate should. Where the two corpora disagree by 5x
# (text_figure_refs born-digital: 34.0% vs 6.6% — plate-based monographs
# against modern articles) the ceiling necessarily sits near the looser
# one, and a corpus like viburnum should tighten it locally via
# CORPUS_SOFT_RATE_CEILINGS rather than run against a ceiling built for
# somebody else's material.
#
# references_match_corpus_papers is the outlier and is deliberately set
# where it can only catch catastrophe. Its premise — that a paper cites
# other papers *in this corpus* — measures corpus cohesion, not reference
# parsing, and cohesion is not a pipeline property: viburnum was assembled
# by mining references outward, so a paper matching nothing in-corpus is
# its normal case, at 52.7% / 84.4%. Total parser failure would still show
# as ~100%, which 90/95 catches; anything subtler this check cannot see on
# a loosely-coupled corpus and should not pretend to.
#
# That reasoning is now confirmed rather than argued (PLAN v1.2 §3). The
# gold transcriptions carry the reference sections as printed, so reference
# *parsing* can be measured directly: of 659 references Grobid extracted
# across the gold corpuscle, **94.8% have a title that appears in the gold
# text of that document**. Grobid is not inventing references, and the
# residue falls on the same axes as everything else —
#
#     Ahuja 2026, Mańko 2020, Hosia 2024, Stepanjants 2014   100%
#     Totton 1965a                                            95.5%
#     Benasso 1976                                            83.9%
#     Eschscholtz 1825 (Fraktur)                              50.0%
#     Linnaeus 1735                                            0.0%
#
# So the loose ceiling above is justified: parsing is good, and the thing
# this check actually varies with is how tightly the corpus cites itself.
# Tightening it would catch corpus assembly, not a pipeline regression.
#
# The OCR rows were measured before the #176 fix and should tighten once
# affected papers are re-ingested.
_SOFT_RATE_CEILINGS = {
    "title_appears_in_text":                {"born_digital": 0.30, "ocr": 0.45},
    "text_figure_refs_have_extracted_figures": {"born_digital": 0.50, "ocr": 0.80},
    "first_author_surname_in_text":         {"born_digital": 0.25, "ocr": 0.25},
    "references_match_corpus_papers":       {"born_digital": 0.90, "ocr": 0.95},
    "alphabet_ratio":                       {"born_digital": 0.12, "ocr": 0.20},
    "extracted_figures_referenced_in_text": {"born_digital": 0.25, "ocr": 0.25},
    "no_duplicate_chunks":                  {"born_digital": 0.05, "ocr": 0.05},
}

def _load_ceilings():
    """Ceilings, with a per-corpus override from ``CORPUS_SOFT_RATE_CEILINGS``.

    The defaults above are calibrated across two corpora, but rates are a
    property of the *material* — how many papers are offprints without title
    pages, how many are plate-based monographs, how tightly the corpus cites
    itself — so a corpus of a different shape needs its own numbers, and
    #185 asked for exactly that. A corpus that sits well inside the defaults
    should tighten them here rather than leave itself covered by a ceiling
    sized for somebody else's material. Point the env var at a JSON file::

        {"title_appears_in_text": {"born_digital": 0.40, "ocr": 0.70}}

    Keys not mentioned keep their default, so an override file only has to
    carry the rows that differ.
    """
    override = os.environ.get("CORPUS_SOFT_RATE_CEILINGS")
    if not override:
        return _SOFT_RATE_CEILINGS
    data, err = _load_json_safe(Path(override))
    if err:
        raise RuntimeError(f"CORPUS_SOFT_RATE_CEILINGS={override}: {err}")
    merged = {k: dict(v) for k, v in _SOFT_RATE_CEILINGS.items()}
    for check, buckets in (data or {}).items():
        if check not in merged:
            raise RuntimeError(
                f"CORPUS_SOFT_RATE_CEILINGS: unknown check {check!r}; "
                f"known: {', '.join(sorted(merged))}"
            )
        merged[check].update(buckets)
    return merged


_FILE_TYPE_MAP = None


def _file_type_bucket(doc_hash):
    """`born_digital` or `ocr` for a document.

    `scanned` and `broken_text_layer` share a bucket: both are OCR'd, both
    show the same premise drift, and there are only 14 of the latter — too
    few for a rate of its own to mean anything.
    """
    global _FILE_TYPE_MAP
    if _FILE_TYPE_MAP is None:
        _FILE_TYPE_MAP = {}
        for h, d in _DOCS:
            data, _ = _load_json_safe(d / "scan_detection.json")
            ft = (data or {}).get("file_type")
            _FILE_TYPE_MAP[h] = "born_digital" if ft == "born_digital" else "ocr"
    return _FILE_TYPE_MAP.get(doc_hash, "ocr")


def _assert_soft_rate(check_name, predicate):
    """Run ``predicate`` over every document and assert the failure rate.

    ``predicate(doc_hash, doc_dir)`` returns True (consistent), False
    (inconsistent) or None (not applicable — missing input, too little text,
    a known-untestable shape). Rates are computed per file-type bucket over
    the applicable documents only.
    """
    if len(_DOCS) < _SOFT_RATE_MIN_DOCS:
        pytest.skip(
            f"corpus is {len(_DOCS)} document(s); a rate needs at least "
            f"{_SOFT_RATE_MIN_DOCS} to mean anything"
        )
    ceilings = _load_ceilings()[check_name]
    buckets = {k: {"n": 0, "failed": []} for k in ceilings}
    for doc_hash, doc_dir in _DOCS:
        try:
            verdict = predicate(doc_hash, doc_dir)
        except Exception as e:                       # noqa: BLE001
            pytest.fail(f"{check_name} raised on {doc_hash}: {e!r}")
        if verdict is None:
            continue
        b = buckets[_file_type_bucket(doc_hash)]
        b["n"] += 1
        if not verdict:
            b["failed"].append(doc_hash)

    problems = []
    report = []
    for name, ceiling in ceilings.items():
        n = buckets[name]["n"]
        failed = buckets[name]["failed"]
        if not n:
            continue
        rate = len(failed) / n
        report.append(f"{name}: {len(failed)}/{n} = {rate:.1%} (ceiling {ceiling:.0%})")
        if rate > ceiling:
            problems.append(
                f"{name}: {rate:.1%} of {n} applicable documents failed "
                f"{check_name}, above the {ceiling:.0%} ceiling. "
                f"First offenders: {', '.join(failed[:10])}"
            )
    assert not problems, "\n".join(problems) + "\n\nAll buckets: " + "; ".join(report)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(params=_DOC_IDS)
def doc(request):
    """Yields (doc_hash, doc_dir)."""
    return request.param, _doc_dir(request)


@pytest.fixture(params=_DOC_IDS)
def summary(request):
    d = _doc_dir(request)
    data, err = _load_json_safe(d / "summary.json")
    assert err is None, f"{request.param}: summary.json parse error: {err}"
    return request.param, data


@pytest.fixture(params=_DOC_IDS)
def metadata(request):
    d = _doc_dir(request)
    return request.param, _load_json_or_skip(d, "metadata.json"), d


@pytest.fixture(params=_DOC_IDS)
def figures_and_chunks(request):
    """Yields (doc_hash, figures_data, chunks_data, doc_dir)."""
    d = _doc_dir(request)
    figs = _load_json_or_skip(d, "figures.json")
    chunks = _load_json_or_skip(d, "chunks.json")
    return request.param, figs, chunks, d


@pytest.fixture(params=_DOC_IDS)
def text_and_metadata(request):
    """Yields (doc_hash, text_data, metadata_data, scan_detection_data)."""
    d = _doc_dir(request)
    text = _load_json_or_skip(d, "text.json")
    md = _load_json_or_skip(d, "metadata.json")
    sd_path = d / "scan_detection.json"
    sd = None
    if sd_path.exists():
        sd, _ = _load_json_safe(sd_path)
    return request.param, text, md, sd


@pytest.fixture(params=_DOC_IDS)
def references_data(request):
    d = _doc_dir(request)
    return request.param, _load_json_or_skip(d, "references.json")


@pytest.fixture(params=_DOC_IDS)
def chunks_data(request):
    d = _doc_dir(request)
    return request.param, _load_json_or_skip(d, "chunks.json")


# ===================================================================
# STRUCTURAL CHECKS
# ===================================================================

EXPECTED_FILES = [
    "summary.json",
    "metadata.json",
    "text.json",
    "figures.json",
    "references.json",
    "chunks.json",
    "scan_detection.json",
]


class TestFilesExist:
    """Every document directory should contain the core JSON files."""

    @pytest.mark.parametrize("doc_hash", _DOC_IDS)
    @pytest.mark.parametrize("filename", EXPECTED_FILES)
    def test_expected_file_exists(self, doc_hash, filename):
        path = _DOC_MAP[doc_hash] / filename
        assert path.exists(), f"{doc_hash}/{filename} missing"


class TestJsonParseable:
    """All core JSON files should parse without error."""

    @pytest.mark.parametrize("doc_hash", _DOC_IDS)
    @pytest.mark.parametrize("filename", EXPECTED_FILES)
    def test_json_parses(self, doc_hash, filename):
        path = _DOC_MAP[doc_hash] / filename
        if not path.exists():
            pytest.skip(f"{filename} missing")
        _, err = _load_json_safe(path)
        assert err is None, f"{doc_hash}/{filename}: {err}"


class TestSummary:
    def test_status_is_success(self, summary):
        doc_hash, data = summary
        status = data.get("processing_summary", {}).get("status")
        assert status == "success", f"{doc_hash}: status={status}"

    def test_no_processing_errors(self, summary):
        doc_hash, data = summary
        errors = data.get("processing_summary", {}).get("errors", [])
        assert len(errors) == 0, (
            f"{doc_hash}: {len(errors)} errors: {errors[:3]}"
        )

    def test_has_pdf_hash(self, summary):
        doc_hash, data = summary
        assert data.get("pdf_hash") == doc_hash

    def test_has_processing_steps(self, summary):
        doc_hash, data = summary
        steps = data.get("processing_summary", {}).get("processing_steps", [])
        assert len(steps) > 0, f"{doc_hash}: no processing steps recorded"


# ===================================================================
# FIGURE ↔ TEXT CROSS-REFERENCING
# ===================================================================

class TestFigureTextConsistency:
    """Check bidirectional consistency between extracted figures and text."""


    def test_figure_chunk_crossrefs_are_symmetric(self, figures_and_chunks):
        """referenced_in_chunks in figures.json and figure_refs in chunks.json
        should agree.

        If a figure claims it's referenced in chunk_5, then chunk_5's
        figure_refs should list that figure's ID, and vice versa.
        """
        doc_hash, figs_data, chunks_data, _ = figures_and_chunks
        figs = figs_data.get("figures", [])
        chunks = chunks_data.get("chunks", [])

        # Build both directions
        fig_claims_chunks = {}  # figure_id -> set of chunk_ids
        for f in figs:
            fid = f.get("figure_id")
            if fid:
                fig_claims_chunks[fid] = set(f.get("referenced_in_chunks", []))

        chunk_claims_figs = {}  # chunk_id -> set of figure_ids
        for ch in chunks:
            cid = ch.get("chunk_id")
            if cid:
                chunk_claims_figs[cid] = set(ch.get("figure_refs", []))

        asymmetric = []
        # Check figure -> chunk direction
        for fid, claimed_chunks in fig_claims_chunks.items():
            for cid in claimed_chunks:
                if cid in chunk_claims_figs and fid not in chunk_claims_figs[cid]:
                    asymmetric.append(f"fig {fid} claims {cid}, but {cid} doesn't list {fid}")

        # Check chunk -> figure direction
        for cid, claimed_figs in chunk_claims_figs.items():
            for fid in claimed_figs:
                if fid in fig_claims_chunks and cid not in fig_claims_chunks[fid]:
                    asymmetric.append(f"chunk {cid} claims {fid}, but {fid} doesn't list {cid}")

        assert not asymmetric, (
            f"{doc_hash}: {len(asymmetric)} asymmetric cross-refs: "
            f"{asymmetric[:5]}"
        )

    def test_figure_files_exist_on_disk(self, figures_and_chunks):
        """Every figure in figures.json should have its PNG on disk."""
        doc_hash, figs_data, _, doc_dir = figures_and_chunks
        figs = figs_data.get("figures", [])
        missing = []
        for fig in figs:
            fname = fig.get("filename")
            if not fname:
                continue
            if not (doc_dir / "figures" / fname).exists():
                missing.append(fname)
        assert not missing, (
            f"{doc_hash}: {len(missing)} figure files missing: {missing[:5]}"
        )


# ===================================================================
# CITATION GRAPH CONSISTENCY
# ===================================================================

class TestCitationGraph:
    """Cross-reference each paper's reference list against the corpus.

    These tests build a corpus-wide index on first use (adds a few seconds).
    They verify that reference extraction is producing metadata consistent
    enough to match known papers.
    """


    def test_paper_is_cited_by_others(self, metadata):
        """Papers published before 2015 in this focused corpus should be cited
        by at least one other paper.

        A paper that exists in the corpus but is never referenced by any other
        corpus paper may have incorrect metadata (wrong year, garbled author
        name) making it unmatchable.
        """
        doc_hash, md, _ = metadata
        year = md.get("year")
        authors = md.get("authors", [])
        if year is None or not authors:
            pytest.skip("missing year or authors")
        if year >= 2015:
            pytest.skip("too recent — may not be cited yet")

        first = authors[0]
        surname = first.get("surname", "") if isinstance(first, dict) else str(first)
        surname = surname.strip().lower()
        if not surname:
            pytest.skip("empty first-author surname")

        # Search all other papers' reference lists for this (surname, year)
        cited = False
        for other_hash, other_dir in _DOCS:
            if other_hash == doc_hash:
                continue
            refs_path = other_dir / "references.json"
            if not refs_path.exists():
                continue
            refs_data, _ = _load_json_safe(refs_path)
            if refs_data is None:
                continue
            for ref in refs_data.get("references", []):
                ref_year = ref.get("year")
                ref_authors = ref.get("authors", [])
                if ref_year != year or not ref_authors:
                    continue
                ref_first = ref_authors[0] if isinstance(ref_authors[0], str) else ""
                ref_surname = ref_first.strip().split()[-1].lower() if ref_first.strip() else ""
                if ref_surname == surname:
                    cited = True
                    break
            if cited:
                break

        # This is informational — many legitimate papers won't match due to
        # OCR errors, transliteration differences, etc. We only flag it
        # as a warning via a soft assertion.
        if not cited:
            pytest.xfail(
                f"{doc_hash}: {surname.title()} ({year}) not cited by any "
                f"other corpus paper — may indicate metadata mismatch"
            )


# ===================================================================
# METADATA PLAUSIBILITY
# ===================================================================

class TestMetadataPlausibility:
    """Cross-check metadata fields against each other and the text."""

    def test_has_title(self, metadata):
        doc_hash, data, _ = metadata
        title = data.get("title", "")
        assert title and title.strip(), f"{doc_hash}: empty title"

    def test_has_year(self, metadata):
        doc_hash, data, _ = metadata
        assert data.get("year") is not None, f"{doc_hash}: missing year"

    def test_year_is_plausible(self, metadata):
        doc_hash, data, _ = metadata
        year = data.get("year")
        if year is None:
            pytest.skip("no year")
        assert 1700 <= year <= 2026, f"{doc_hash}: implausible year {year}"

    def test_has_authors(self, metadata):
        doc_hash, data, _ = metadata
        assert len(data.get("authors", [])) > 0, f"{doc_hash}: no authors"

    def test_filename_year_matches_metadata_year(self, metadata):
        """When the filename contains a 4-digit year, it should match the
        extracted metadata year. Disagreement suggests extraction error."""
        doc_hash, data, _ = metadata
        filename = data.get("filename", "")
        md_year = data.get("year")
        if md_year is None:
            pytest.skip("no metadata year")
        m = re.search(r"(\d{4})", filename)
        if not m:
            pytest.skip("no year in filename")
        fn_year = int(m.group(1))
        if fn_year < 1700 or fn_year > 2026:
            pytest.skip("filename number is not a plausible year")
        assert fn_year == md_year, (
            f"{doc_hash}: filename year {fn_year} != metadata year {md_year} "
            f"(file: {filename})"
        )


# ===================================================================
# TEXT QUALITY SIGNALS
# ===================================================================

class TestTextQuality:
    """Detect text extraction problems via statistical signals."""

    def test_has_text(self, text_and_metadata):
        doc_hash, text_data, _, _ = text_and_metadata
        body = text_data.get("text", "")
        assert len(body) > 0, f"{doc_hash}: empty text"

    def test_text_min_length(self, text_and_metadata):
        doc_hash, text_data, _, _ = text_and_metadata
        body = text_data.get("text", "")
        assert len(body) >= 100, (
            f"{doc_hash}: suspiciously short text ({len(body)} chars)"
        )

    def test_chars_per_page(self, text_and_metadata):
        """Papers with very few characters per page likely have extraction
        failures (e.g. OCR missed pages, or docling returned empty cells)."""
        doc_hash, text_data, _, _ = text_and_metadata
        body = text_data.get("text", "")
        pages = text_data.get("pages")
        if not pages or pages == 0:
            pytest.skip("no page count")
        cpp = len(body) / pages
        # A typical page has 2000-4000 chars. Below 200 is suspicious.
        assert cpp >= 200, (
            f"{doc_hash}: only {cpp:.0f} chars/page across {pages} pages — "
            f"likely extraction failure"
        )


    def test_has_page_count(self, text_and_metadata):
        doc_hash, text_data, _, _ = text_and_metadata
        pages = text_data.get("pages")
        assert pages is not None and pages > 0, (
            f"{doc_hash}: missing or zero page count"
        )


# ===================================================================
# CHUNK QUALITY
# ===================================================================

class TestChunkQuality:
    """Check chunking consistency."""

    def test_has_chunks(self, chunks_data):
        doc_hash, data = chunks_data
        ch = data.get("chunks", [])
        assert len(ch) > 0, f"{doc_hash}: no chunks"

    def test_few_empty_chunks(self, chunks_data):
        doc_hash, data = chunks_data
        ch = data.get("chunks", [])
        if not ch:
            pytest.skip("no chunks")
        empty = [i for i, c in enumerate(ch) if not c.get("text", "").strip()]
        assert len(empty) <= max(1, len(ch) * 0.1), (
            f"{doc_hash}: {len(empty)}/{len(ch)} chunks have empty text"
        )


    def test_chunk_length_not_degenerate(self, chunks_data):
        """Flag papers where most chunks are extremely short (< 50 chars) —
        suggests the chunker is splitting on every line."""
        doc_hash, data = chunks_data
        ch = data.get("chunks", [])
        if len(ch) < 3:
            pytest.skip("too few chunks")
        lengths = [len(c.get("text", "")) for c in ch]
        short = sum(1 for l in lengths if l < 50)
        frac = short / len(lengths)
        assert frac <= 0.5, (
            f"{doc_hash}: {frac:.0%} of chunks are < 50 chars — "
            f"chunker may be over-splitting"
        )


# ===================================================================
# SOFT CONSISTENCY CHECKS — asserted as corpus rates (#185)
# ===================================================================
#
# Each predicate is the per-paper logic that used to be an assertion,
# returning True / False / None (not applicable). The per-paper detail is
# not lost: a rate that breaches its ceiling names the first offenders, so
# triage still starts from a document hash.


def _check_extracted_figures_referenced_in_text(doc_hash, doc_dir):
    """Numbered figures should be mentioned somewhere in the text."""
    figs_data, _ = _load_json_safe(doc_dir / "figures.json")
    chunks_data, _ = _load_json_safe(doc_dir / "chunks.json")
    if figs_data is None or chunks_data is None:
        return None
    numbered = [
        f for f in figs_data.get("figures", [])
        if f.get("figure_number") is not None
        and f.get("figure_type") in ("figure", "plate")
    ]
    if not numbered:
        return None
    all_text = " ".join(ch.get("text", "") for ch in chunks_data.get("chunks", []))
    text_nums = _extract_text_figure_numbers(all_text)
    unreferenced = [
        f["figure_number"] for f in numbered
        if not _fig_number_digit_keys(f["figure_number"]) & text_nums
    ]
    return len(unreferenced) / len(numbered) <= 0.5


def _check_text_figure_refs_have_extracted_figures(doc_hash, doc_dir):
    """Figure numbers cited in text should have entries in figures.json.

    Plate-based works break this premise routinely: plates are bound
    separately, numbered in their own sequence, or live in another volume.
    """
    figs_data, _ = _load_json_safe(doc_dir / "figures.json")
    chunks_data, _ = _load_json_safe(doc_dir / "chunks.json")
    if figs_data is None or chunks_data is None:
        return None
    all_text = " ".join(ch.get("text", "") for ch in chunks_data.get("chunks", []))
    text_nums = _extract_text_figure_numbers(all_text)
    if not text_nums:
        return None
    extracted = set()
    for f in figs_data.get("figures", []):
        fn = f.get("figure_number")
        if fn is not None:
            extracted |= _fig_number_digit_keys(fn)
    return len(text_nums - extracted) / len(text_nums) <= 0.5


def _check_first_author_surname_in_text(doc_hash, doc_dir):
    """The first author's surname should appear in their own paper."""
    text_data, _ = _load_json_safe(doc_dir / "text.json")
    md, _ = _load_json_safe(doc_dir / "metadata.json")
    if text_data is None or md is None:
        return None
    authors = md.get("authors", [])
    if not authors:
        return None
    first = authors[0]
    surname = (first.get("surname", "") if isinstance(first, dict) else str(first)).strip()
    if not surname or len(surname) < 2:
        return None
    body = text_data.get("text", "")
    if _transliteration_likely(surname, body):
        return None
    return _norm(surname) in _norm(body)


def _check_title_appears_in_text(doc_hash, doc_dir):
    """The .bib title should appear in the body text.

    Weakest premise of the set: for this material the recorded title is
    frequently curated, translated or normalised relative to what the page
    prints, and offprints and extracts often have no title page at all.
    """
    text_data, _ = _load_json_safe(doc_dir / "text.json")
    md, _ = _load_json_safe(doc_dir / "metadata.json")
    if text_data is None or md is None:
        return None
    title = md.get("title", "").strip()
    if not title or len(title) < 10:
        return None
    body = text_data.get("text", "")
    if _transliteration_likely(title, body):
        return None
    title_norm, body_norm = _norm(title), _norm(body)
    if title_norm in body_norm:
        return True
    return title_norm[:min(40, len(title_norm))] in body_norm


def _check_alphabet_ratio(doc_hash, doc_dir):
    """Low alphabetic fraction suggests OCR garbage."""
    text_data, _ = _load_json_safe(doc_dir / "text.json")
    if text_data is None:
        return None
    body = text_data.get("text", "")
    if len(body) < 100:
        return None
    return sum(1 for c in body if c.isalpha()) / len(body) >= 0.40


def _check_no_duplicate_chunks(doc_hash, doc_dir):
    """Exact duplicate chunk text within a paper suggests a chunking bug."""
    data, _ = _load_json_safe(doc_dir / "chunks.json")
    if data is None:
        return None
    ch = data.get("chunks", [])
    if len(ch) < 2:
        return None
    texts = [c.get("text", "").strip() for c in ch if c.get("text", "").strip()]
    seen, dupes = set(), 0
    for t in texts:
        if t in seen:
            dupes += 1
        seen.add(t)
    return dupes <= max(1, len(texts) * 0.05)


def _check_references_match_corpus_papers(doc_hash, doc_dir):
    """A paper with many references should match at least one corpus paper."""
    data, _ = _load_json_safe(doc_dir / "references.json")
    if data is None:
        return None
    refs = data.get("references", [])
    if len(refs) < 15:
        return None
    ref_authors = " ".join(
        (r.get("authors") or [""])[0] if isinstance((r.get("authors") or [""])[0], str)
        else ""
        for r in refs
    )
    if ref_authors.strip() and _substantially_nonlatin(ref_authors):
        return None
    index = _build_corpus_index()
    for ref in refs:
        year, authors = ref.get("year"), ref.get("authors", [])
        if not authors or year is None:
            continue
        first = authors[0] if isinstance(authors[0], str) else ""
        surname = first.strip().split()[-1].lower() if first.strip() else ""
        if surname and any(h != doc_hash for h in index.get((surname, year), ())):
            return True
    return False


class TestSoftConsistencyRates:
    """Corpus-health rates for checks whose per-paper premise is unreliable.

    One assertion per check instead of one per paper. See the rationale and
    the ceiling table above `_SOFT_RATE_CEILINGS`.
    """

    def test_title_appears_in_text_rate(self):
        _assert_soft_rate("title_appears_in_text", _check_title_appears_in_text)

    def test_text_figure_refs_have_extracted_figures_rate(self):
        _assert_soft_rate("text_figure_refs_have_extracted_figures",
                          _check_text_figure_refs_have_extracted_figures)

    def test_first_author_surname_in_text_rate(self):
        _assert_soft_rate("first_author_surname_in_text",
                          _check_first_author_surname_in_text)

    def test_references_match_corpus_papers_rate(self):
        _assert_soft_rate("references_match_corpus_papers",
                          _check_references_match_corpus_papers)

    def test_alphabet_ratio_rate(self):
        _assert_soft_rate("alphabet_ratio", _check_alphabet_ratio)

    def test_extracted_figures_referenced_in_text_rate(self):
        _assert_soft_rate("extracted_figures_referenced_in_text",
                          _check_extracted_figures_referenced_in_text)

    def test_no_duplicate_chunks_rate(self):
        _assert_soft_rate("no_duplicate_chunks", _check_no_duplicate_chunks)

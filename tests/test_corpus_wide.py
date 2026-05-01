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

# ---------------------------------------------------------------------------
# Discovery and loading
# ---------------------------------------------------------------------------

_DEFAULT_OUTPUT = Path("/nfs/roberts/project/pi_cwd7/cwd7/output")


def _output_dir():
    override = os.environ.get("CORPUS_OUTPUT_DIR")
    if override:
        return Path(override)
    repo_output = Path(__file__).parent.parent / "output"
    if repo_output.is_dir():
        return repo_output
    return _DEFAULT_OUTPUT


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

    def test_extracted_figures_referenced_in_text(self, figures_and_chunks):
        """Figures in figures.json should be mentioned somewhere in the text.

        Figures with a figure_number that never appear in any chunk's text
        suggest either a spurious figure extraction or a text extraction gap.
        Graphical elements and unclassified figures are excluded — they often
        lack explicit text references.
        """
        doc_hash, figs_data, chunks_data, _ = figures_and_chunks
        figs = figs_data.get("figures", [])
        chunks = chunks_data.get("chunks", [])

        # Numbered figures/plates only (skip graphical_element, unclassified)
        numbered = [
            f for f in figs
            if f.get("figure_number") is not None
            and f.get("figure_type") in ("figure", "plate")
        ]
        if not numbered:
            pytest.skip("no numbered figures/plates")

        # Collect all figure numbers mentioned in chunk text
        all_text = " ".join(ch.get("text", "") for ch in chunks)
        text_nums = _extract_text_figure_numbers(all_text)

        unreferenced = [
            f["figure_number"] for f in numbered
            if not _fig_number_digit_keys(f["figure_number"]) & text_nums
        ]
        # Allow some slack — captions can use unusual formatting
        frac = len(unreferenced) / len(numbered) if numbered else 0
        assert frac <= 0.5, (
            f"{doc_hash}: {len(unreferenced)}/{len(numbered)} numbered figures "
            f"not referenced in text: {unreferenced[:10]}"
        )

    def test_text_figure_refs_have_extracted_figures(self, figures_and_chunks):
        """Figure numbers cited in running text should have corresponding
        entries in figures.json.

        If the text says 'Fig. 5' but no figure 5 was extracted, the figure
        extraction pipeline may have missed it.
        """
        doc_hash, figs_data, chunks_data, _ = figures_and_chunks
        figs = figs_data.get("figures", [])
        chunks = chunks_data.get("chunks", [])

        all_text = " ".join(ch.get("text", "") for ch in chunks)
        text_nums = _extract_text_figure_numbers(all_text)
        if not text_nums:
            pytest.skip("no figure references in text")

        extracted_keys = set()
        for f in figs:
            fn = f.get("figure_number")
            if fn is not None:
                extracted_keys |= _fig_number_digit_keys(fn)

        missing = text_nums - extracted_keys
        # Many papers reference figures in other papers, so be lenient.
        # But if >50% of referenced figures are missing, something is wrong.
        frac = len(missing) / len(text_nums) if text_nums else 0
        assert frac <= 0.5, (
            f"{doc_hash}: {len(missing)}/{len(text_nums)} figure numbers cited "
            f"in text but not extracted: {sorted(missing)[:10]}"
        )

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

    def test_references_match_corpus_papers(self, references_data):
        """Some fraction of each paper's references should match other papers
        in the corpus (by first-author surname + year).

        This is a soft check — many references will be to papers outside this
        siphonophore corpus. But a paper with 30 references and zero corpus
        matches is suspicious (likely broken reference parsing).
        """
        doc_hash, data = references_data
        refs = data.get("references", [])
        if len(refs) < 3:
            pytest.skip("too few references to check")

        index = _build_corpus_index()
        matches = 0
        for ref in refs:
            year = ref.get("year")
            authors = ref.get("authors", [])
            if not authors or year is None:
                continue
            # First author surname — refs store authors as plain strings
            first_author = authors[0] if isinstance(authors[0], str) else ""
            # Extract surname (last whitespace-separated token)
            surname = first_author.strip().split()[-1].lower() if first_author.strip() else ""
            if surname and (surname, year) in index:
                # Don't count self-matches
                corpus_matches = [h for h in index[(surname, year)] if h != doc_hash]
                if corpus_matches:
                    matches += 1

        # We don't assert a minimum match rate because some papers cite
        # entirely outside the corpus. Instead, we collect the data — the
        # test passing with 0 matches is fine, but the output is useful
        # for identifying reference parsing issues at scale.
        # Flag only if there are many refs and literally zero matches.
        if len(refs) >= 15:
            assert matches >= 1, (
                f"{doc_hash}: 0/{len(refs)} references match any corpus paper "
                f"(by first-author + year). Likely broken reference parsing."
            )

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

    def test_first_author_surname_in_text(self, text_and_metadata):
        """The first author's surname should appear somewhere in the paper text.

        Authors virtually always appear in their own paper (title page, running
        header, or self-citations). Absence suggests garbled metadata.
        """
        doc_hash, text_data, md, _ = text_and_metadata
        authors = md.get("authors", [])
        if not authors:
            pytest.skip("no authors")
        first = authors[0]
        surname = first.get("surname", "") if isinstance(first, dict) else str(first)
        surname = surname.strip()
        if not surname or len(surname) < 2:
            pytest.skip("first author surname too short or empty")
        body = text_data.get("text", "")
        assert surname.lower() in body.lower(), (
            f"{doc_hash}: first author surname '{surname}' not found in text"
        )

    def test_title_appears_in_text(self, text_and_metadata):
        """The paper title (or a substantial substring) should appear in the
        extracted text. Missing title suggests either wrong title in metadata
        or a text extraction gap on the first page."""
        doc_hash, text_data, md, _ = text_and_metadata
        title = md.get("title", "").strip()
        if not title or len(title) < 10:
            pytest.skip("title too short to check")
        body = text_data.get("text", "")
        # Try full title first, then first 40 chars (titles can get truncated)
        title_lower = title.lower()
        body_lower = body.lower()
        if title_lower in body_lower:
            return  # pass
        # Try a substantial prefix
        prefix = title_lower[:min(40, len(title_lower))]
        assert prefix in body_lower, (
            f"{doc_hash}: title not found in text. "
            f"Title: '{title[:60]}...'"
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

    def test_alphabet_ratio(self, text_and_metadata):
        """Text with a low fraction of alphabetic characters is likely OCR
        garbage or extraction artifacts."""
        doc_hash, text_data, _, _ = text_and_metadata
        body = text_data.get("text", "")
        if len(body) < 100:
            pytest.skip("text too short")
        alpha = sum(1 for c in body if c.isalpha())
        ratio = alpha / len(body)
        # Normal text is ~70-80% alphabetic. Below 40% is suspicious.
        assert ratio >= 0.40, (
            f"{doc_hash}: only {ratio:.1%} alphabetic characters — "
            f"likely OCR noise or extraction artifacts"
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

    def test_no_duplicate_chunks(self, chunks_data):
        """Exact duplicate chunk text within a paper suggests a chunking bug."""
        doc_hash, data = chunks_data
        ch = data.get("chunks", [])
        if len(ch) < 2:
            pytest.skip("too few chunks")
        texts = [c.get("text", "").strip() for c in ch if c.get("text", "").strip()]
        seen = set()
        dupes = []
        for t in texts:
            if t in seen:
                dupes.append(t[:60])
            seen.add(t)
        assert len(dupes) <= max(1, len(texts) * 0.05), (
            f"{doc_hash}: {len(dupes)} duplicate chunks: {dupes[:3]}"
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

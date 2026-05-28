# Testing & Evaluation

Two complementary test suites check pipeline output quality. Both run against
already-processed output — they never re-run the pipeline.

1. **Ground-truth tests** — per-paper checks against human-curated answers
   (3 papers currently). High precision, narrow coverage.
2. **Corpus-wide tests** — structural and content consistency checks across
   all 1,787 papers. No ground truth needed — they verify things that can be
   checked programmatically.

## Quick start

```bash
conda activate corpus

# All tests
python -m pytest tests/ -v

# Ground-truth tests only
python -m pytest tests/test_text_extraction.py tests/test_figure_extraction.py tests/test_reference_extraction.py -v

# Corpus-wide tests only
python -m pytest tests/test_corpus_wide.py -v

# Corpus-wide, one category
python -m pytest tests/test_corpus_wide.py -v -k "TestFigureTextConsistency"

# Compact failure output (useful for corpus-wide)
python -m pytest tests/test_corpus_wide.py -v --tb=line
```

## What's tested

### Ground-truth tests (per-paper)

| Module | What it checks |
|--------|---------------|
| `test_text_extraction.py` | Required phrases present, garbage absent, minimum text length |
| `test_figure_extraction.py` | Figure count, specific figures (number/caption/type), files on disk |
| `test_reference_extraction.py` | Reference count, specific references, metadata (title/year/authors) |

### Corpus-wide content consistency tests

| Test class | What it checks |
|------------|---------------|
| `TestFilesExist` | All 7 core JSON files present in every document directory |
| `TestJsonParseable` | All core JSON files parse without error |
| `TestSummary` | Processing status is "success", no errors recorded, hash matches |
| `TestFigureTextConsistency` | Bidirectional figure ↔ text cross-referencing (see below) |
| `TestCitationGraph` | Reference lists match other corpus papers; pre-2015 papers are cited (see below) |
| `TestMetadataPlausibility` | Year/author/title cross-checks against text and filename (see below) |
| `TestTextQuality` | Chars/page ratio, alphabet fraction, minimum text length |
| `TestChunkQuality` | Duplicate chunks, empty chunks, over-splitting detection |

#### Figure ↔ text cross-referencing

Three checks for bidirectional consistency:

- **Extracted figures referenced in text**: numbered figures/plates in
  `figures.json` should be mentioned in chunk text (e.g. "Fig. 1").
  Unreferenced figures may be spurious extractions or graphical noise.
- **Text figure refs have extracted figures**: figure numbers cited in running
  text (parsed via regex for "Fig.", "Figure", "Plate" patterns) should have
  corresponding entries in `figures.json`. Missing figures suggest extraction
  gaps.
- **Cross-ref symmetry**: `referenced_in_chunks` in figures.json and
  `figure_refs` in chunks.json must agree in both directions. Asymmetry is a
  pipeline bug.

#### Citation graph consistency

- **References match corpus papers**: for each paper's reference list, how
  many references match another corpus paper by (first-author-surname, year)?
  Papers with 15+ references and zero corpus matches likely have broken
  reference parsing.
- **Paper is cited by others**: papers published before 2015 in this focused
  siphonophore corpus should be cited by at least one other corpus paper.
  Uncited papers may have garbled metadata making them unmatchable. Uses
  `pytest.xfail` (expected failure) since legitimate misses exist.

#### Metadata plausibility

- **Filename year vs metadata year**: when the PDF filename contains a 4-digit
  year, it must match the GROBID-extracted year. Disagreement flags extraction
  error.
- **First author surname in text**: the first author's surname should appear
  somewhere in the paper (title page, headers, or self-citations). Absence
  suggests garbled metadata.
- **Title appears in text**: the full title (or a 40-char prefix) should appear
  in the extracted text. Missing title suggests wrong metadata or first-page
  extraction gap.

#### Text quality signals

- **Chars per page**: papers below 200 chars/page likely have OCR or extraction
  failures on some pages.
- **Alphabet ratio**: text with < 40% alphabetic characters is likely OCR
  garbage or extraction artifacts.

## Running subsets

```bash
# One paper
python -m pytest tests/ -v -k "Dunn"

# One quality dimension
python -m pytest tests/test_figure_extraction.py -v

# Against a different output directory
CORPUS_OUTPUT_DIR=/path/to/output python -m pytest tests/ -v
```

## Ground truth files

Ground truth lives in `tests/ground_truth/*.yaml`, one file per paper. Tests
auto-discover all YAML files in this directory — no Python changes needed to
add papers or assertions.

### Format

```yaml
pdf_filename: Example.pdf
hash: af043530e5dd          # 12-char SHA-256 prefix

metadata:
  title: "Paper Title"       # checked as substring (case-insensitive)
  year: 2005                 # exact match
  authors_contain:           # each surname must appear
    - surname: Smith
    - surname: Jones

text:
  must_contain:              # phrases that must appear in extracted text
    - "key taxonomic term"
    - "important anatomical structure"
  must_not_contain:          # OCR artifacts or garbage strings
    - "ÿÿÿ"
  min_length: 5000           # minimum character count

figures:
  min_count: 3               # at least this many figures extracted
  expected:                  # specific figures to verify
    - figure_number: 1
      caption_contains: "species name"
      figure_type: figure    # figure | plate | graphical_element | unclassified

references:
  min_count: 10              # at least this many references extracted
  expected:                  # specific references to verify
    - title_contains: "keyword"
      authors_contain: ["Smith"]
      year: 2001
```

Any section or field can be omitted — absent sections are skipped, not failed.

### Adding a new assertion

Edit a YAML file and add a line. For example, to test that Dunn 2005 mentions
"bract" in the text:

```yaml
text:
  must_contain:
    - "Marrus claudanielis"
    - "bract"               # ← add this line
```

### Adding a new paper

1. Create `tests/ground_truth/<filename_stem>.yaml`
2. Set `hash` to the paper's 12-char SHA-256 prefix (find it in the output
   directory name, or run `sha256sum demo/Paper.pdf | cut -c1-12`)
3. Fill in whatever ground truth you know
4. Run `python -m pytest tests/ -v -k "<filename_stem>"`

## Current ground truth papers

| Paper | Hash | Type | Tests |
|-------|------|------|-------|
| Dunn-etal2005_Marrus | `af043530e5dd` | Modern, English | text, figures, references, metadata |
| Siebert_etal2011 | `45d2af65e152` | Modern, English | text, figures, references, metadata |
| Pugh2001_Erenna | `4fe914163f59` | Plates, English | text, figures, references, metadata |

### Candidates to add

| Paper | Hash | Challenge | Expected value |
|-------|------|-----------|---------------|
| Schneider1891 | `ef8482d9cb44` | German Fraktur, scanned | OCR quality regression detection |
| Quoy_Gaimard1827 | `01091787348c` | Historical French | Multilingual OCR |
| Alekseev1984 | `b756815902e7` | Russian, scanned | Cyrillic OCR |
| Stepanjants1970 | `dde93d15a5e8` | Broken text layer | Scan detection, forced OCR |
| Pages_etal1991 | `3eafb0775ece` | Modern English | Baseline coverage |

## Design notes

- Tests read from the main output directory by default
  (`/nfs/roberts/project/pi_cwd7/cwd7/output`). Override with
  `CORPUS_OUTPUT_DIR`.
- The conftest computes SHA-256 hashes from `demo/*.pdf` to map filename
  stems to output directories. The `hash` field in YAML overrides this.
- `pytest.skip()` is used when output files are missing or ground truth
  sections are absent, so partial pipeline runs don't cause false failures.

# Testing & Evaluation

Two complementary test suites check pipeline output quality. Both run against
already-processed output — they never re-run the pipeline.

1. **Ground-truth tests** — per-paper checks against human-curated answers.
   High precision, narrow coverage.
2. **Corpus-wide tests** — structural and content consistency checks across
   the current corpus. No ground truth needed — they verify things that can be
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

### Scoped figure downloads and reverse proxies

`tests/test_signed_figure_urls.py` covers scope tampering, expiry, restart,
licensing rechecks, and isolation from MCP bearer authentication. The opt-in
live test runs the route from `deploy/nginx.conf` in an unprivileged, disposable
Linux Docker container and fetches whole figures and panels over real HTTP:

```bash
docker pull nginx:stable-alpine
CORPUS_TEST_NGINX=1 pytest -q tests/test_signed_figure_urls.py
```

It binds only high loopback ports, removes its container, and does not modify
an existing deployment. This complements, rather than replaces, full MCP smoke
testing against a read-only real bundle. Record the image digest with release
evidence when running this acceptance test.

### Embedding update and recovery tests

`tests/test_embedding_updates.py` uses real temporary LanceDB tables and a
deterministic two-dimensional backend, so it needs no GPU or model download.
Run it with `tests/test_package_for_serve.py` to check document replacement,
content/metadata invalidation, legacy duplicate repair, empty-document pruning,
interrupted writes, no-op resume, whole-index bundle validation and equality of
clean/incremental logical rows. These run in T0. They verify storage/update
semantics, not embedding quality or full upstream BibTeX/config invalidation;
the latter still needs the end-to-end gate in the
[update contract](OVERVIEW.md#corpuscle-update-contract).

`tests/test_metadata_resume.py` drives both Stage 1 resume gates with real
BibTeX parsing, metadata extraction and completion records while stubbing the
expensive PDF/figure work. It checks per-paper edits, entry addition/removal,
renames with and without a bib match, identical-copy path updates, dry-run
non-mutation, preservation of existing vision results, and metadata/chunk/
reference equality with a clean build.

`tests/test_config_updates.py` exercises the same two gates for configured
consumer/descendant changes and default restoration. It compares rechunking
and panel-disable updates with clean builds; checks reset of renamed/split
images, model changes, raster refresh and stale Docling sidecars; verifies TEI
input/payload validation and archival when Grobid is absent; and injects
publication failures, stage interruption and standalone vision failures.
It also checks CPU/vision handoff preservation, read-only configuration drift
and CLI precedence. PDF/model work is stubbed, so these are deterministic
update/recovery tests, **not caption-quality scores or full corpus acceptance**.
External-service/model provenance and broader clean/incremental parity remain
separate release gates.

`tests/test_grobid_recovery.py` uses a controllable fake service with real
metadata parsing and both resume gates. It covers disabled/enabled transitions,
startup outages, individual request and parse failures, malformed service
responses, preserved valid caches, version/declared-producer changes, unknown
versions, strict-network errors, URL precedence, applied timeouts and network-free
dry runs. It compares recovered metadata/reference artifacts against a clean
build and checks that status reads persisted outcomes after unrelated reruns.
It does not test a live Grobid deployment or verify custom model contents.

### Page-level visual audit

When a score or acceptance prompt points to a particular page, generate the
optional report from that document's build artifacts:

```bash
python -m pipeline.page_report /path/to/output/documents/<HASH> \
  --pages 12-13,17
```

`page_report.html` is one self-contained, locally viewable file. It shows the
rendered `processed.pdf` page beside selectable Docling text and provides
toggleable overlays for PDF word cells, extracted figure boxes, projected ROI
boxes, and chosen/rejected caption candidates. Page statistics and the full
caption evidence trail make coordinate or ownership errors inspectable without
manually joining the PDF, `docling_doc.json`, `figures.json`, and text output.
Normal builds do not generate it, and served-bundle distillation excludes it.
The default safety limit is 200 selected pages; select a smaller range or pass
`--max-pages 0` deliberately for a longer document.

### Reference reconciliation audit

After rebuilding the bibliographic authority database, audit the missing-work
ranking against its raw observation evidence:

```bash
python tools/qc/reference_reconciliation.py \
  --db /path/to/output/biblio_authority.sqlite \
  --min-citations 2 --limit 50 --out reference-reconciliation.json
```

The command is read-only and requires the v1.3 observation schema. It reports
the current, mapped and unmapped observation/citing-document populations,
compatibility citation edges, mapping methods and producer versions, then gives
bounded raw/parsed evidence for each ranked missing work. A substantive
same-title/year corpus candidate is a review signal only; `author_set_match`
identifies the narrower exact identity evidence used by the deterministic
resolver. The report never merges works and does not change the MCP response
surface.

### Caption-binding and panel-split fidelity

Score a built corpuscle against the independent page transcriptions:

```bash
python tools/qc/caption_binding.py \
  --gold /path/to/transcriptions \
  --corpuscle /path/to/output \
  --out caption-binding.json
```

The read-only report separates the raw figure record from the default MCP
evidence types. Figure-number recall/precision measures same-page ownership;
`reported_pair_capacity_rate` gives the best recall possible if each page keeps
its current count of distinct reported number pairs and those pairs are
relabelled perfectly. This diagnoses missing upstream number evidence; it is
not a theoretical ceiling on a rebuild that discovers additional labels. The
panel section separately reports declaration recall/precision, exact
letter sets, and label recall/precision. Both gold parsers are independent of
the production functions they measure, and numeric figures sharing a plate
are excluded from the letter-panel denominator. Declaration precision uses
same-page, same-number gold figure blocks without panel enumerations as its
negative set. Page diagnostics retain the
expected, reported, missing and surplus numbers and panel labels needed to
open the corresponding page report.

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
| `TestCitationGraph` | Reference lists match other corpus papers; reference-evaluation fixtures check cited-paper coverage (see below) |
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
- **Paper is cited by others**: the reference evaluation checks that older
  in-scope papers are cited by at least one other corpus paper. This is a
  collection-specific consistency expectation, not a requirement for every
  corpus. Uncited papers may have garbled metadata making them unmatchable;
  the check uses `pytest.xfail` because legitimate misses exist.

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

prompts:                     # citation-grounding round-trips (#79)
  - id: marrus_synopsis
    prompt: "Summarize the discovery of Marrus claudanielis with
             citations from the corpus."
    expected_citations:      # each must be emitted via format_citations
      - work_id: "corpus:..."
      - work_id: "10.1234/..."
    forbidden_hallucinations:  # negative trip-wires
      - author_surname: "Schneider"
        attributed_journal_not: "Bull. Mus. Comp. Zool."
    expected_failure: false  # set true for known gaps; xfails
```

Any section or field can be omitted — absent sections are skipped, not failed.

### Citation-grounding tests (`prompts:` block)

Each entry in the `prompts:` block runs as a real Claude API round-trip
with the `format_citations` MCP tool wired in as a callable. The
harness (`tests/test_prompt_quality.py`) asserts that every expected
`work_id` is emitted via the tool (the `formatted` output must appear
verbatim in the response) and that no parenthetical-citation tokens
appear that don't trace to a logged tool call — the latter catches the
recombination class that #79 exists to prevent.

`forbidden_hallucinations` is the negative trip-wire: each entry's
`author_surname` and `attributed_journal_not` must not co-occur within
a ~200-character window of the response. Use it to pin known incidents
(the original taxonomist-feedback case being the seed) as a
no-regression guard.

**Gating.** These tests cost real API calls and are skipped unless
`RUN_PROMPT_QUALITY=1`. Within CI, they're intended for the
release-time lane with `ANTHROPIC_API_KEY` in the runner secrets —
not on every PR. Locally:

```bash
RUN_PROMPT_QUALITY=1 ANTHROPIC_API_KEY=sk-... \
  python -m pytest tests/test_prompt_quality.py -v
# Optional: override the model (default claude-haiku-4-5-20251001)
PROMPT_QUALITY_MODEL=claude-sonnet-4-6-20251001 ...
```

The harness's scoring functions (`_score`, regex matchers) have
in-process unit tests that always run regardless of the gate, so a
regression in the harness itself can't ride into a release.

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

## Build regression references

Use the operator-side snapshot alongside independent gold scoring:

```bash
python -m tools.qc.build_reference snapshot /path/to/baseline --out baseline.json
python -m tools.qc.build_reference snapshot /path/to/candidate --out candidate.json
python -m tools.qc.build_reference compare baseline.json candidate.json --out comparison.json
```

Snapshot and comparison outputs are created exclusively: an existing file is
an error, never silently replaced. Comparison exits 1 for semantic differences,
document additions/removals, missing primary artifacts or hard failures; this
means **review required**, not necessarily a regression. Unchanged warning
populations are listed separately and still need an explicit disposition.
There is no automatic baseline acceptance or blanket warning suppression.

The snapshot hashes primary JSON content, retaining schema versions and
reporting producer versions separately. Build-root paths are relocated and
top-level producer stamps do not count as content changes. Reference fields
are retained by ordinal so an unchanged count cannot conceal a damaged title,
author list, DOI or raw citation. A reordered bibliography also needs review.
Ground-truth reference expectations may require `raw_contains` as well as
`title_contains` and `authors_contain`: source preservation and correct field
parsing are separate assertions. Gold-only fixtures with explicit document
hashes are skipped when that document is absent from the demo build.
Prepared-PDF and TEI byte hashes are diagnostic: generated IDs and PDF metadata
can differ without an extraction change. They are not standalone failure gates.

Schema v2 also fingerprints bibliography tables (including current observation
mappings, membership, permissions and curation provenance), taxon mentions,
taxonomy, and all logical LanceDB rows. Vector reads stream in batches and
compare exact float values and row multiplicity per document; physical row
order and transaction generation IDs do not matter. Decoded RGBA image hashes
catch pixel changes while ignoring PNG compression/metadata and legacy query
crops. Database bookkeeping timestamps and mention row IDs are excluded;
whether a work was curator-imported is retained. Historical raw observations
and reconciliation decisions are append-only history, not current-build
equality targets. Current authority-table differences still require review,
even when caused by retained historical knowledge.

Missing optional databases/indexes are explicit `null`, not a proof of their
completeness; check that the expected indexes exist for the acceptance corpus.
Schema v1 references must be regenerated from both retained builds, never
silently upgraded. Agreement with a baseline does not establish correctness:
keep the independent gold reports and review their known misses. Malformed
JSON or SQLite is a hard error, not an omitted paper.

## Design notes

- Tests read from the main output directory by default
  (`/nfs/roberts/project/pi_cwd7/cwd7/output`). Override with
  `CORPUS_OUTPUT_DIR`.
- The conftest computes SHA-256 hashes from `demo/*.pdf` to map filename
  stems to output directories. The `hash` field in YAML overrides this.
- `pytest.skip()` is used when output files are missing or ground truth
  sections are absent, so partial pipeline runs don't cause false failures.

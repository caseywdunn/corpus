# PDF QC workflow

How to find bad papers in a corpus, decide whether to keep them, and
reprocess after fixing root causes. Backs
[#54](https://github.com/caseywdunn/corpus/issues/54).

## The two gates

`corpus` distinguishes two reasons a paper might be excluded from the
served bundle:

- **Skip flag** — the operator decides this paper isn't worth serving.
  Curated via the BibTeX `serve = {false}` field, persisted to
  `works.serve` in `biblio_authority.sqlite`. Honored at deploy time
  by `mcpsrv.bundle` (`corpus run`'s bundle distillation step). The
  per-paper artifacts stay in the build bundle so the paper still
  appears in `corpus status` rollups.
- **Publishable gate** — the figure-licensing decision (#51). See
  [LICENSING.md](LICENSING.md). Orthogonal to the skip flag: the skip
  gate decides "appears in MCP at all"; the publishable gate decides
  "image is reusable in derived publications".

This file is about the skip flag.

## The triage workflow

### 1. Find the worst papers

```bash
corpus status --sort-by quality_flag_count --tail 20
corpus status --sort-by stage_failure_count --tail 20
```

Lists the 20 papers carrying the most quality_flag entries (or stage
failures). The `quality_flag` set comes from v0.2's silent-failure
detectors (#36): `empty_text`, `low_text_density`,
`gibberish_after_ocr`, `ocr_pages_blanked`, `ocr_no_text_recovered`,
`zero_references_unexpected`, `single_token_chunks`, `all_black_figures`.
`empty_text` ignores Docling's `<!-- image -->` layout placeholders;
`ocr_no_text_recovered` is an error when OCR exits successfully but leaves
every output page without a text layer. Per-paper QC metrics
(citation_resolution_rate, dictionary_hit_rate, headfoot_pollution_score,
…) are a #54 follow-up; gate counts are the proxy until those land.

### 2. Decide which to skip

```bash
corpus status --propose-skips --min-flags 2
```

Lists papers carrying ≥2 quality flags that aren't yet `serve = false`,
with the gate names that fired and a paste-ready BibTeX block. The
output looks like:

```bibtex
% af043530e5dd  3 flags: empty_text, gibberish_after_ocr, zero_references_unexpected
serve = {false},
```

Paste the `serve = {false},` line (and optionally a `servereason`)
into the matching BibTeX entry, then:

```bash
corpus bib import my_corpus.bib
```

`works.serve` flips to 0; the next `corpus run` excludes the paper
from the served bundle.

### 3. Curate the reason

The `servereason` BibTeX field carries a short tag from a closed
vocabulary (free-text fallback OK):

| reason | meaning |
|---|---|
| `content-no-value` | bibliographic-only, pure tables, conference program |
| `quality-irrecoverable` | OCR irrecoverable, scan damaged, language pack absent |
| `out-of-scope` | wrong topic, ended up in the corpus by mistake |
| `duplicate` | superseded by another entry |
| `rights-restricted` | not publishable; ties into [#51](https://github.com/caseywdunn/corpus/issues/51) |

```bibtex
@article{smith2020,
  ...
  serve        = {false},
  servereason  = {content-no-value: 40-page taxonomic key, no novel data},
}
```

### 4. Compare against precedent

```bash
corpus status --skipped
```

Groups currently-excluded papers by `serve_reason` so you can see what
the team has already decided is out of scope. Use this when triaging
a new candidate — does it look like what you've already skipped?

### 5. Fix root causes + reprocess

When a quality flag fires because of an *operator-fixable* problem
(missing language pack, broken Grobid, etc.), fix the root cause then:

```bash
corpus run --re-process-flagged gibberish_after_ocr
```

This clears `pipeline_state.json` for every paper carrying the named
gate, so the per-stage resume logic re-extracts them on the upcoming
run. Other papers stay cached. The `--re-process-flagged` flag is
typed as a single gate per invocation; chain runs for multiple gates.

### `ocr_pages_blanked` is the one that needs a config change first

Most gates describe the document. This one describes the *machine that
built it*, so re-running without changing anything reproduces it — on a
different, load-dependent set of documents. Comparing two full builds of
the same library, 31 documents lost >80% of their text while 28
independently gained it back (#254).

The cause is worker oversubscription: ocrmypdf sizes its Tesseract pool
from `multiprocessing.cpu_count()`, which reports the whole machine even
inside an 8-CPU allocation. Every page then runs on a sliver of a core,
crosses `ocr.tesseract_page_timeout`, and is copied into the output blank
while ocrmypdf exits 0. The pipeline now derives `--jobs` from the
allocation rather than the host, so on a build from v1.2 onward this gate
firing means either an explicit `ocr.jobs` that is too high or a
genuinely slow page. On an older build, or where the allocation cannot be
read, set it by hand to match the CPUs the job requested:

```yaml
ocr:
  jobs: 8          # match --cpus-per-task
```

then `corpus run --re-process-flagged ocr_pages_blanked`. Raising
`ocr.tesseract_page_timeout` treats the symptom — the page is still
running at 1/8 speed, just for longer.

## Skipped papers in the system

| Layer | What happens |
|---|---|
| Per-paper artifacts | Unaffected. summary.json, figures.json, etc. stay on disk |
| `corpus status` rollups | Skipped papers still appear (so you can compare against new candidates) |
| Build bundle | Skipped papers stay |
| Served bundle | Skipped papers excluded by `mcpsrv.bundle` (#54) |
| MCP `get_papers` / `list_papers` | Skipped papers don't appear (reads served bundle) |
| MCP citation graph | A skipped paper can still appear as a *citation target* of a non-skipped paper (authority record exists) |

This is the documented contract. Reverting a skip is non-destructive:
flip the BibTeX field, re-run `corpus bib import`, re-run
`corpus run` (the bundle step picks it back up).

## Granularity

Skip is **paper-level only** in v0.3. No per-figure skip mechanism
exists — defer until a real case (one bad figure in an otherwise fine
paper) forces it.

## What's missing in v0.3

- New per-paper QC metrics (citation_resolution_rate,
  dictionary_hit_rate, page_text_presence_histogram,
  figure_to_page_ratio, caption_density, headfoot_pollution_score,
  grobid_header_confidence). Adds richer signal for `--sort-by`.
  Tracked as polish; gate counts cover the worst-first triage today.
- Lexicon-coverage section in `corpus status --report` (per #57's
  reference script). Defer to a polish pass.

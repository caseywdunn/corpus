# PLAN.md — Corpus pipeline (v0.2)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing →
MCP-serving stack on a ~1,800-paper siphonophore corpus. v0.2 hardens
v0.1's rough edges and fills in the deferred design items: vision pass
at corpus scale, geographic mention layer, batch-script cleanup, and
the granular-resume restructure now under evaluation.

This document is scoped to v0.2 work. Architectural background and
pipeline internals live in [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md);
per-feature history in [CHANGELOG.md](CHANGELOG.md); HPC operations in
[dev_docs/BOUCHET.md](dev_docs/BOUCHET.md); deployment in
[dev_docs/DEPLOY.md](dev_docs/DEPLOY.md). Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues); the
v0.2-scoped subset is listed below.

## 1. v0.2 punch list

### Hardening — close out v0.1 rough edges

- [#9](https://github.com/caseywdunn/corpus/issues/9) — Install
  `deu_latf` Fraktur pack on Bouchet. ~6–8 19th-c. German scans
  (Goldfuss 1820, Pagenstecher 1869, Brandt 1837, Dönitz 1871, …)
  currently OCR to whitespace; reprocess them after the pack lands.
- [#16](https://github.com/caseywdunn/corpus/issues/16) — Figure-number
  extraction on old/scanned papers. ~538/1787 papers have `Fig. N`
  mentions in body text but no extracted figure number. Improve
  `parse_figure_number` heuristics or fall back to running-text
  mention positions when caption parsing fails.
- [#17](https://github.com/caseywdunn/corpus/issues/17) — MCP server
  identity. Partially addressed in v0.2-dev: `__version__` from
  [version.py](version.py) now surfaces via the `bundle_info` tool.
  Remaining: server name should be corpuscle-aware (e.g.,
  `corpus:siphonophores` rather than `corpus`).
- [#6](https://github.com/caseywdunn/corpus/issues/6) — `deploy/stack.yaml`
  default-VPC assumption. Fails on accounts without a default VPC;
  parameterize or pre-flight.
- [#22](https://github.com/caseywdunn/corpus/issues/22) —
  `ingest_taxonomy.py` DwC field-name fix (`taxon` vs. `taxa`).
- [#25](https://github.com/caseywdunn/corpus/issues/25) — `pngquant`
  missing from `environment.yaml`.
- Bouchet batch-script cleanup
  ([#19](https://github.com/caseywdunn/corpus/issues/19),
  [#20](https://github.com/caseywdunn/corpus/issues/20),
  [#21](https://github.com/caseywdunn/corpus/issues/21)) — old code in
  `batch_embed.sh`, library-load problems across batch scripts, CUDA
  / torch fixes.
- [#7 part 3](https://github.com/caseywdunn/corpus/issues/7) — Any
  remaining work on the in-text citation graph beyond parts 1 and 2
  already merged.

### Vision + figures

- [#11](https://github.com/caseywdunn/corpus/issues/11) — **Vision pass
  at corpus scale.** Pipeline is wired (`batch_pass3b.sh`,
  Qwen2.5-VL-7B local backend, Claude remote backend) and tested on
  the demo set; v0.1 didn't run it on the full corpus. Resolves the
  bulk of the 6,841 pre-existing `missing_figures[]` records and the
  compound-figure splits (Pass 3c).
- [#27](https://github.com/caseywdunn/corpus/issues/27) — Add SLURM
  array processing to `batch_pass3b.sh` so it parallelizes the same
  way `batch_pipeline.sh` does. Prerequisite for #11 at any reasonable
  wallclock.

### Indices + features

- [#13](https://github.com/caseywdunn/corpus/issues/13) —
  **Geographic mention layer.** Largest new design item; see [§2](#2-geographic-mention-layer-issue-13).
- [#5](https://github.com/caseywdunn/corpus/issues/5) — **Streamable
  HTTP transport with OAuth**, replacing the SSE-with-bearer-token
  transport that v0.1 ships. Streamable HTTP is the current MCP
  standard; bearer-token auth is fine for ~20 collaborators but not
  for wider distribution.
- [#26](https://github.com/caseywdunn/corpus/issues/26) — Round-trip
  bibliography: MCP export + shell import.
- [#23](https://github.com/caseywdunn/corpus/issues/23) — Expand
  taxonomy ingest beyond marine taxa (currently WoRMS-shaped
  assumptions); plant DwC archives as the test case.
- [#24](https://github.com/caseywdunn/corpus/issues/24) — Expand
  anatomy lexicon.

### Internal

- [#15](https://github.com/caseywdunn/corpus/issues/15) — Refactor
  `mcp_server.py`. The 2050-line single file has outgrown a single
  module; split per-concern (papers / taxonomy / bibliography /
  figures / transports). No behavior change.
- [#12](https://github.com/caseywdunn/corpus/issues/12) — Stage 2
  parallelism. `batch_embed.sh` already parallelizes; assess whether
  further work is needed (per-document chunking concurrency, BGE-M3
  batch sizing).

## 2. Geographic mention layer (issue #13)

Bibliographic and taxonomic mention layers shipped in v0.1; geographic
is the third and last layer in the per-corpuscle `*_mentions` family.
Same shape as the others — external table that points into chunk text,
rebuilt independently of the pipeline (see
[dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) for the existing layers).

### Schema

```sql
CREATE TABLE geo_mentions (
    mention_id     INTEGER PRIMARY KEY,
    locality_id    TEXT,                   -- FK to localities (NULL if unresolved)
    corpus_hash    TEXT NOT NULL,
    chunk_index    INTEGER NOT NULL,
    char_start     INTEGER NOT NULL,
    char_end       INTEGER NOT NULL,
    mention_text   TEXT NOT NULL,
    latitude       REAL,
    longitude      REAL,
    depth_m        REAL,
    confidence     REAL DEFAULT 1.0,
    method         TEXT NOT NULL           -- 'spacy_ner' | 'geonames_match' | 'regex_coords'
);

CREATE TABLE localities (
    locality_id    TEXT PRIMARY KEY,
    name           TEXT NOT NULL,
    latitude       REAL,
    longitude      REAL,
    country        TEXT,
    geonames_id    INTEGER,
    source         TEXT NOT NULL           -- 'geonames' | 'manual' | 'inferred'
);
```

### Extraction

Sampling localities appear in this literature in several forms: named
places ("Villefranche-sur-Mer", "Bay of Naples"), coordinates
("32°15′N, 64°40′W"), station references ("Station 42", "Sta. 1247"),
depth ranges ("200–400 m"), and vessel names ("R/V *Dana*",
"HMS *Challenger*").

1. **Named-location NER** — spaCy `GPE`/`LOC` entities as a starting
   point. Filter for geographic relevance (not all GPEs are sampling
   locations).
2. **Geocoding** — resolve named locations to coordinates via the
   GeoNames API or a local GeoNames dump. Record `geonames_id` for
   traceability.
3. **Coordinate extraction** — regex for degree-minute-second and
   decimal-degree patterns. Associate with nearby taxon mentions for
   taxon × locality co-occurrence.
4. **Structured sampling events (deferred).** Extract (location,
   depth, date, vessel, taxon) tuples from methods sections. Harder;
   may benefit from LLM extraction. Out of scope for v0.2.

### MCP tool surface

| Tool | Purpose |
| --- | --- |
| `get_locality_mentions(locality, radius_km=None)` | All text spans mentioning a locality (or within radius) |
| `get_taxon_localities(name)` | Co-occurring taxon × locality pairs across the corpus |
| `get_locality_taxa(locality)` | All taxa mentioned near a given locality |

### Feasibility

Highest-effort and lowest-initial-precision of the three mention
layers. Named locations are detectable but associating them with the
right taxon — and distinguishing sampling locations from other
geographic references ("the Mediterranean fauna") — requires context.
Start with named locations + coordinates only; structured sampling
events are a later enrichment.

## 3. Proposal under evaluation: granular per-stage resume (#28)

[#28](https://github.com/caseywdunn/corpus/issues/28) proposes
replacing the all-or-nothing `summary.json` completion marker with
per-stage status, plus renaming `process_corpus.py` and the numbered
passes to describe-what-they-do (`extract.py`, `parse_metadata.py`,
…). Motivation: selective redo (re-run only Grobid metadata after a
Grobid outage; re-run only annotation after a taxonomy bump) without
re-paying for Docling extraction. Decision needed before this lands
in v0.2 — restructure work touches every batch script, the resume
logic, and the MCP server's startup checks.

## 4. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`. v0.1 made Q1, Q2, Q4–Q8 answerable; Q3 is the only
one still gated on unbuilt indices (trait extraction, deferred to
post-v0.2 — [#14](https://github.com/caseywdunn/corpus/issues/14)).

| # | Pattern | Status after v0.1 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer (§2) for the locality side |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred (#14, post-v0.2) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Needs vision pass at corpus scale (#11) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Needs vision-pass figure tagging at corpus scale (#11) |

## 5. Versioning + release ritual

`__version__` in [version.py](version.py) is the single source of
truth and is stamped into every persistent artifact (bundle manifest,
MCP `bundle_info`). [CONTRIBUTING.md](CONTRIBUTING.md) covers the
branching model and release ritual; the short version: `dev` carries
a PEP 440 pre-release suffix (`0.2.0.dev0` → `0.2.0a1` → `0.2.0`),
the release commit drops the suffix, the next commit on `dev`
reintroduces one for the next target.

## 6. Out of scope for v0.2

- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to
  warrant its own plan section when picked up; v0.3 candidate.
- Multi-region failover, autoscaling, or blue/green deploys for the
  AWS served bundle. Single-instance is fine until it isn't.
- Authentication beyond bearer tokens / OAuth (Cognito, institutional
  SSO).
- Mirror to Cloudflare R2 / Backblaze B2 for cost. Defer until S3
  egress shows up on a bill.
- A thin HTML/web UI on top of the MCP server. Out of scope until the
  MCP-only experience has actual non-Claude-Desktop users.

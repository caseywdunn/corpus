# PLAN.md — Corpus pipeline (v1.0)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing
→ MCP-serving stack. v0.2 hardened internals (per-stage resume, input
fingerprints, structured failure schema, quality gates). v0.3 collapsed
the user surface into one CLI plus a per-corpuscle `config.yaml`. v0.4
closed the silent-failure modes external testers and clean-room HPC
hosts surfaced, and gated the next cycle on tiered CI. v0.5 was the
served-bundle quality cycle — LLM-side citation amalgamation
([#79](https://github.com/caseywdunn/corpus/issues/79)), the dossier
suite ([#76](https://github.com/caseywdunn/corpus/issues/76)), and a
pinned ML stack after an unpinned docling release silently broke macOS
arm64 extraction ([#98](https://github.com/caseywdunn/corpus/issues/98),
[#99](https://github.com/caseywdunn/corpus/issues/99)).

**v0.6 (2026-06-04) was the API-freeze cycle and it landed.** The MCP
surface is frozen at **38 tools** with a uniform error payload and a
consistent pagination convention: the redundant singular tools were
removed in favor of batched plurals, `lexicon_matrix` went
summary-by-default, `get_citation_graph` gained breadth caps, the
output-type profile ontology
([#101](https://github.com/caseywdunn/corpus/issues/101)) shipped, and
`tests/test_freeze_contract.py` now guards the surface mechanically.
Ops work landed alongside: `/healthz` capability reporting
([#91](https://github.com/caseywdunn/corpus/issues/91)), per-invocation
run logs ([#90](https://github.com/caseywdunn/corpus/issues/90)),
`corpus debug-pdf` ([#92](https://github.com/caseywdunn/corpus/issues/92)),
and the `--figure-panels` selector
([#102](https://github.com/caseywdunn/corpus/issues/102)). The v0.6
punch list is preserved in the
[v0.6.0 tag's history](https://github.com/caseywdunn/corpus/blob/v0.6.0/dev_docs/PLAN.md).

**v1.0 is the installability cycle.** The surface is settled; what is
not settled is whether a stranger can install this software and have it
work. Two external install reports say not reliably —
[#145](https://github.com/caseywdunn/corpus/issues/145) (macOS Apple
Silicon, blocked at Grobid after a ~20 GB download) and
[#153](https://github.com/caseywdunn/corpus/issues/153) (HPC, hand-built
cache redirection and a non-localhost Grobid URL) — and as of this
writing **`conda env create -f environment.yaml` produces a broken
environment** ([#156](https://github.com/caseywdunn/corpus/issues/156)).

The organizing principle: **a green CI badge must mean that a fresh
install works** — and, just as important, that a *silent* badge doesn't
mean anything at all.

> **Correction to #156 and #158, established 2026-07-30.** Both issues
> attribute the invisible drift to a cached conda env. That is wrong:
> `setup-miniconda` does not cache the solved environment and this repo
> adds no `actions/cache` for it, so T0/T1/T2 re-resolve
> `environment.yaml` on every run. The only caches are integration.yml's
> HuggingFace model-hub caches. The real mechanism is **timeliness**: CI
> runs only on push, and `dev` sat untouched from 2026-07-06 to
> 2026-07-30. `mcp` 2.0.0 shipped in that window and broke every clean
> install for 24 days; the first push after it (`f6cff9a`) failed
> immediately with the expected 18 collection errors. This makes pinning
> *more* important than #158 argued — with a fresh solve every run, an
> unpinned dep makes CI itself nondeterministic, so any push can break
> for reasons unrelated to the push — and it makes the fix cheaper: a
> `schedule:` trigger, not a cache-key redesign.

The gaps are therefore a missing clock (now
[`clean-room.yml`](../.github/workflows/clean-room.yml), T3) and a
documented install path no job exercises — nothing ever runs
`docker-compose.yml` ([#161](https://github.com/caseywdunn/corpus/issues/161)).
Every other item in this cycle is downstream: until a clean install is
verified continuously, no fix in it can be trusted to survive.

Secondary through-line: **no silent wrongness.** The 1.0 bug list is
dominated by code that produces plausible-looking bad output rather than
errors — a bib parser that discards 88% of a file with one warning
([#141](https://github.com/caseywdunn/corpus/issues/141)), curated
references stuck in the warning tier
([#142](https://github.com/caseywdunn/corpus/issues/142) →
[#152](https://github.com/caseywdunn/corpus/issues/152)), taxonomy
skipped into empty `taxa.json`
([#139](https://github.com/caseywdunn/corpus/issues/139)), OCR falling
back to English on Cyrillic pages
([#160](https://github.com/caseywdunn/corpus/issues/160)). These are the
failures a new user cannot debug, and 1.0 is the version new users get.

**No new MCP tools.** The v0.6 freeze holds. Where 1.0 changes served
behavior it is to make an existing tool honest
([#154](https://github.com/caseywdunn/corpus/issues/154),
[#143](https://github.com/caseywdunn/corpus/issues/143)), not to widen
the surface.

Writing the formal 1.0 API-stability *policy* doc — deliberately
deferred out of v0.6 — belongs to this cycle's release step.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](../DEPLOY.md);
platform-portability criteria in
[PLATFORM_SMOKE.md](PLATFORM_SMOKE.md). Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## 0. Release gate

One gate, and it is not a checklist item — it is the definition of the
cycle:

> **A clean-room install from `environment.yaml` must be verified by CI,
> not by hand, before 1.0 is tagged.**

Wave 0 built it:
[`.github/workflows/clean-room.yml`](../.github/workflows/clean-room.yml)
(**T3**) runs conda env → `pip install -e .` → the real
`docker-compose.yml` → demo `corpus run` → bundle → `corpus_required` →
SSE round-trip, with the HuggingFace cache off so a first-run model
download is exercised for real.

**How to actually invoke it.** GitHub honors `schedule:` and
`workflow_dispatch:` only for workflows present on the **default branch**
(`main`), so while T3 lives on a feature branch it is unregistered —
`gh workflow run clean-room.yml` returns `HTTP 404` and no weekly run
fires. A gate that cannot be invoked until after the release merge is not
a gate, so the lane also triggers on:

* **a pull request targeting `main`** — which *is* the release proposal,
  so the gate runs automatically on it rather than depending on someone
  remembering to dispatch. This is the intended path: it satisfies the
  gate *before* the merge, and dovetails with CONTRIBUTING.md's existing
  "wait for CI on main to be green before tagging" step.
* **a push that modifies `clean-room.yml` itself** — so a change to the
  lane is validated by running it, from any branch.

Once 1.0 is on `main`, the weekly `schedule:` takes over as the standing
drift detector.

[`dev_docs/ec2_smoke.sh`](ec2_smoke.sh) (**T3-bare**) stays manual and
pre-release: it covers the one thing T3 cannot, the bare-host bootstrap
(apt, miniforge install) on a real Ubuntu EC2 instance, against
[PLATFORM_SMOKE.md](PLATFORM_SMOKE.md)'s criteria.

Everything else in §1 is verified *against* T3.

No open branches or pre-merge gates. (`issue-corpus-run-hpc-1-engine`,
gated in the v0.6 plan on compute-node acceptance tests, is merged into
`dev`; the remote branch can be deleted.)

## 1. v1.0 punch list

Waves are ordered by dependency, not by severity. Wave 0 must land
first because nothing else can be verified on a clean install until it
does. Within a wave, items are independent and parallelizable across
issue branches.

### Wave 0 — restore installability

- [x] **Migrate `mcpsrv/` to the `mcp` 2.x API and pin exactly**
  ([#156](https://github.com/caseywdunn/corpus/issues/156)). `mcp` is
  unpinned and now resolves to `2.0.0`, which removed
  `mcp.server.fastmcp`. Two import sites break —
  `mcpsrv/app.py:31` (`FastMCP`) and `mcpsrv/tools/figures.py:24`
  (`Image`) — yielding 18 test-collection errors and a `corpus serve`
  that cannot start. **Decision: bump to 2.x and pin `mcp==2.0.0`
  exactly**, not `mcp<2`. The issue body carries a migration inventory
  verified against 2.0.0 by live introspection (~8 lines across 4 files
  plus one test fix). Acceptance **must** include
  `corpus serve --transport sse` end-to-end, since DEPLOY.md's AWS path
  rides SSE and static inspection cannot derisk it.
- [x] **Make dependency drift visible**
  ([#158](https://github.com/caseywdunn/corpus/issues/158)). Two parts,
  both required:
  (a) **A scheduled clean-room lane** —
  [`.github/workflows/clean-room.yml`](../.github/workflows/clean-room.yml),
  weekly + `workflow_dispatch`, so a clean-install failure surfaces in
  days rather than whenever someone next pushes. This is the §0 gate.
  Per the correction above it is a *clock*, not a cache fix; it also
  disables the HuggingFace cache so a genuine first-run model download
  is exercised, and drives the real `docker-compose.yml`. The bare-host
  half of the old T3 (apt + miniforge bootstrap on EC2) stays manual as
  **T3-bare**.
  (b) **Pin the remaining nine pip deps** with a one-line rationale
  each, as #98 did for the ML stack: `ocrmypdf`, `lancedb`,
  `langdetect`, `mcp`, `uvicorn`, `anthropic`, `pytesseract`,
  `qwen-vl-utils`, `accelerate`. The issue recommends pinning only the
  two proven-volatile ones; at 1.0 the cost of a bump ritual is lower
  than the cost of a stranger's broken install, so pin all nine. Both
  `environment.yaml` and `requirements.txt`, per CONTRIBUTING.md
  §"Dependencies — two files, on purpose".

### Wave 1 — first-run experience

Everything a user hits in their first hour. Mostly one-liners and docs;
the exception is #159, which is the largest new surface in the cycle.

- [x] **Fix the Grobid healthcheck**
  ([#157](https://github.com/caseywdunn/corpus/issues/157)).
  `docker-compose.yml:46` probes with `curl`, absent from the
  `lfoppiano/grobid:0.8.1` image, so `corpus-grobid` reports
  `(unhealthy)` forever while working fine. First thing a new user sees
  after `docker compose up -d grobid`. Use the JVM/wget path available
  in the image, or drop to a TCP probe.
- [x] **Exercise `docker-compose.yml` in CI**
  ([#161](https://github.com/caseywdunn/corpus/issues/161)). T1 starts
  Grobid as a GHA `services:` container with different `JAVA_OPTS`, so
  the compose file's image, heap, container name, and healthcheck are
  covered by nothing — which is how #146 and #157 both shipped. Add a
  short job that runs the real file, waits for `/api/isalive`, asserts
  `docker inspect` reports `healthy` (this is what guards #157 from
  regressing), and runs `corpus check`. **Resolved:** the `JAVA_OPTS`
  divergence closed toward CI — `docker-compose.yml` now also passes
  `-XX:-UseContainerSupport`. The T1-compose job caught the NPE on its
  first run, so #72 is not GHA-specific and the compose default was
  broken for users on affected cgroup-v2 hosts. `-Xmx4g` pins the heap
  regardless, so disabling container-aware sizing costs nothing.
- [x] **Drop `do_picture_classification`**
  ([#140](https://github.com/caseywdunn/corpus/issues/140)).
  `pipeline/extract.py:78` sets it `True`, so docling downloads and runs
  `DocumentFigureClassifier-v2.5` on every PDF — but `figure_type` comes
  from our own `classify_figure()` heuristic
  (`pipeline/figures.py:425`), so the output is never read. One line,
  removes an entire HuggingFace download from first run, and it was the
  proximate cause of a CI outage when anonymous runner IPs got 429'd.
  Confirm nothing consumes docling's picture-class annotations before
  flipping, per the issue. Re-runs the figure/extract stage on existing
  corpuscles — changelog it.
- [x] **`corpus prefetch` + model pre-flight**
  ([#159](https://github.com/caseywdunn/corpus/issues/159)). First run
  must reach HuggingFace for docling's layout model, TableFormer, and
  BGE-M3. There is no way to fetch them ahead of time, no way to check
  whether they are cached, and no documentation of where they land — the
  only prefetch logic in the repo is inline Python in
  `integration.yml`, unavailable to users. Add the command (respecting
  `figures.panel_detection` so it doesn't pull the ~16 GB local VLM
  unless asked), a `corpus check` probe, and an INSTALL.md section on
  `HF_HOME` / `TRANSFORMERS_CACHE` / `HF_HUB_OFFLINE`. Have CI call the
  real command instead of maintaining its own copy. Land **after** #140,
  which removes one of the downloads.
- [x] **`corpus check`: validate the OCR toolchain**
  ([#160](https://github.com/caseywdunn/corpus/issues/160)). The only
  `shutil.which` calls live in `pipeline/scan.py` — checked at use time,
  deep in a run. `_cmd_check` (`pipeline/cli.py:869`) never probes
  `tesseract` / `ocrmypdf` / `ghostscript`, nor whether
  `tools/install_tessdata.sh` was run, though the README calls it
  required. Skip that step and Cyrillic OCRs against the English pack
  with only a buried warning. Fail on missing binaries; warn on missing
  language packs, naming the codes and the fixing invocation.
  `_available_tesseract_langs()` already exists.
- [x] **Document the shared-filesystem / HPC install**
  ([#153](https://github.com/caseywdunn/corpus/issues/153)). A user's
  own report: redirect pip/conda caches and `HF_HOME` into project
  space, temporarily reset `$HOME`, run Grobid under Singularity in a
  SLURM job, and point `grobid.url` at the allocated node rather than
  `localhost`. All true, none of it documented. Folds naturally into
  #159's cache section and BOUCHET.md. **Windows/WSL is explicitly not
  supported** and needs no table row.
- [x] **Move test deps out of runtime dependencies**
  ([#162](https://github.com/caseywdunn/corpus/issues/162)).
  `pytest`, `pyflakes`, `ipykernel` are in `[project].dependencies`
  (`pyproject.toml:43-45`), so every install pulls a test runner and a
  Jupyter kernel. Move to `[project.optional-dependencies].dev`; keep
  them unconditional in `environment.yaml`, which *is* the dev env.
- [x] **Bound `requires-python`**
  ([#163](https://github.com/caseywdunn/corpus/issues/163)).
  `pyproject.toml:9` says `>=3.12` with no upper bound while the env
  pins `python=3.12`, CI tests 3.12 only, and the #98 known-good ML set
  was verified against nothing else. Set `>=3.12,<3.13`, or add a 3.13
  lane and widen deliberately — the current unbounded-and-untested
  state is the only unacceptable one.

### Wave 2 — data integrity

Changes stored artifacts, so it should land before users build
corpuscles they intend to keep.

- [x] **Fix the bib parser's silent truncation**
  ([#141](https://github.com/caseywdunn/corpus/issues/141)). A single
  unbalanced brace discards every entry after it: on a 19,834-entry
  export the parse stopped at ~1.75 MB and imported 2,258 entries, with
  one WARNING and a Summary that looked clean. **Three changes, not
  one:**
  (a) **Escape handling** — neither depth counter honors `\{` / `\}`
  (`bib/parser.py:122-130` entry-level, `:62-71` field-level), while the
  quoted-value scanner right below the second one does (`:76-79`). So
  Grobid OCR garbage like `author = {Des, Ej\{aims}` *creates* the
  imbalance. This is the root cause and is a couple of lines.
  (b) **Recovery** at the next top-level `\n@` boundary instead of
  `break` at `:131-133`, so a malformed entry costs one entry.
  (c) **Reconciliation** — count malformed/skipped entries in the
  Summary and warn loudly when the parsed count diverges from the number
  of top-level `@` markers. This is what turns the failure from
  invisible into obvious, and it also covers `_parse_fields`, which has
  no `depth != 0` check at all and silently drops the remaining fields
  of an entry with *no warning whatsoever*.
- [x] **Stamp `bib_imported_at` on matched-but-unchanged entries**
  ([#142](https://github.com/caseywdunn/corpus/issues/142), completing
  [#100](https://github.com/caseywdunn/corpus/issues/100)).
  `bib/importer.py:391-395` `continue`s on the no-change branch before
  reaching `apply_entry()`, making the #100 stamping logic unreachable
  dead code. Cause (a) of #100's three causes was therefore never
  actually fixed in v0.6 — a full round-trip re-import stamped 20 of
  19,834 matched entries. Since an authority-DB rebuild clears
  `bib_imported_at` and forces a re-import, this defeats curation on
  every rebuild. While in the file, collapse the vestigial
  `return 0 if counters["no_match"] == 0 else 0` at `:509` — both
  branches return 0.
- [x] **Verify #152 closes with #142**
  ([#152](https://github.com/caseywdunn/corpus/issues/152)). "Reference
  warnings are overzealous" is a symptom, not an independent bug: with
  `bib_imported_at` unstamped, `CorpusIndex.provenance()` keeps
  returning `grobid_reconciled` and `format_citations` keeps warning on
  hand-curated works. **Acceptance criterion for both:** an unchanged
  `.bib` re-import, followed by a reconcile, leaves every matched work
  in the silent `bib` tier. If warnings persist after #142, the residue
  is #100 cause (b) — `_merge_phase1_into_ghost`
  (`bib/reconcile.py:313-325`) dropping bib fields — and belongs here
  too.
- [x] **Make a missing taxonomy loud**
  ([#139](https://github.com/caseywdunn/corpus/issues/139) follow-ups 1
  and 2). `pipeline/runner.py:305` gates the taxa/anatomy stage on
  `taxonomy_db is not None or bool(lexicons)`, so an absent
  `taxonomy.sqlite` silently skips taxon extraction — on the first full
  Bouchet run that meant 1763 papers with empty `taxa.json` and an empty
  `taxon_mentions.sqlite`. `corpus check` already warns when
  `source: worms` is configured without a built sqlite
  (`pipeline/cli.py:1002-1019`); the gap is at *run* time. Escalate the
  skip to an ERROR naming the pre-build workflow, and surface it in
  `--dry-run`. The docs half (build on a login node, export to DwC-A) is
  done; README §Taxonomy should still note that `source: worms` requires
  internet.

### Wave 3 — served-surface correctness

Making frozen tools honest. No signatures widen; #154 narrows one
response shape, which is the last moment that is cheap.

- [x] **Fix figure licensing**
  ([#154](https://github.com/caseywdunn/corpus/issues/154)). Three
  defects that must be fixed **together**, because §1 and §3 pull in
  opposite directions and both follow from making the gate — not the
  client — responsible:
  (§1) `get_figure` injects license fields regardless of profile
  (`mcpsrv/tools/figures.py:556`), and `get_figure_url` does the same
  (`:860-868`), so a model just *permitted* to display a figure reads
  `"publishable": false` beside it and self-censors. Docstrings actively
  invite this ("so a publication-bound client can self-filter").
  Suppress the fields under `figure_licensing == "permissive"` or gate
  them behind `include_licensing=True`; keep
  `license`/`license_url`/`attribution`, which captions need in every
  profile; add a line to `default_instructions.md` saying the status is
  informational under `report`.
  (§2) `derive_publishable` (`bib/authority.py:401-430`) returns
  `0, "unknown"` both for asserted `all-rights-reserved` and for "no
  license record." In the served corpus that collapse is total: 55,177
  works (86%) are `publishable=0, license_source=unknown` and **not one
  work is asserted all-rights-reserved**. Distinguish the states
  (`public_domain` / `licensed_open` / `restricted` / `undetermined` /
  `no_record`), say which one caused a refusal, and warn at build time
  on unrecognized license strings instead of silently NULLing them.
  (§3) `get_figure_roi_image` (`figures.py:594`) takes no `profile` and
  never calls the gate, so a client refused by `get_figure_image` under
  `manuscript` can obtain the same pixels — including the whole
  uncropped figure when the ROI has no pixel region. Add the parameter
  and apply `_figure_licensing_refusal` before the crop is written.
  **Also fold in** the two `figure_http.py` inconsistencies found in
  review: the route passes the raw `?profile=` to `resolve_profile`
  without validating it via `get_profile` first — contradicting that
  function's own docstring, so a typo'd profile silently resolves to the
  server default instead of erroring — and the hand-rolled query parser
  never URL-decodes, so a label containing `%20` misses its crop.
  `urllib.parse.parse_qs` fixes the second and duplicate-parameter
  handling at once.
- [x] **Expand lexicon synonyms in figure retrieval**
  ([#143](https://github.com/caseywdunn/corpus/issues/143)).
  Ingestion is synonym-aware — `anatomy.json` records every surface form
  against its `canonical` term — but retrieval throws that away and does
  a literal substring count of the one string the caller passed
  (`mcpsrv/tools/figures.py:270` in `get_figures_for_lexicon_term`,
  `:512` in `get_figure_dossier_for_term`). So querying `wing` misses
  captions saying `ala`, querying `ala` misses `wing`, and the tool
  never checks that the term is a known lexicon entry at all. Resolve
  the query through the lexicon to its canonical entry, match against
  the full synonym set, and report which surface form hit. Additive to
  the response; the parameter list does not change.

### Wave 4 — release

- [ ] **Close the v0.6-era review issues.** #146–#151 are all fixed on
  `dev` (`a30f95e`, `1dbb69a`) and were held open pending a version
  bump; #145 is their umbrella and closes with them. Verify each against
  the shipped code rather than the commit message, then close with a
  note to @rdmpage.
- [x] **Write the 1.0 API-stability policy** — deferred out of v0.6 by
  design. What "frozen" means for the 38 tools, what constitutes a
  breaking change, and the deprecation path. Belongs next to
  `dev_docs/MCP_TOOLS.md`. **Shipped** as
  [API_STABILITY.md](API_STABILITY.md), linked from MCP_TOOLS.md and the
  README. It documents the half of the contract that
  `tests/test_freeze_contract.py` already enforces and is explicit about
  what is *not* covered — on-disk artifacts (separately versioned by
  `ARTIFACT_SCHEMA_VERSION`), the CLI, and extracted content, since a
  stable API does not promise stable data and v1.0's re-OCR change is the
  proof.
- [x] **Settle the distribution channel.** PyPI is **not** the right
  primary channel: pip cannot install tesseract, ghostscript, pandoc,
  tectonic, or Grobid, so `pip install` would yield a package that
  imports and then fails on the first scanned PDF — the exact failure
  class this cycle exists to remove. conda-forge *can* carry the native
  toolchain but a feedstock needs every dependency on conda-forge, and
  the pip-only block (`docling`, `lancedb`, `mcp`, `ocrmypdf`, …) is
  not. **Decision: git + conda stays canonical for 1.0**, stated
  explicitly rather than by omission; version the clone
  (`git clone --branch v1.0.0`) and fix `INSTALL.md:138`, which still
  offers `@v0.3.0`. For citability use **Zenodo** (DOI per release via
  the GitHub integration; `CITATION.cff` already exists) — that is what
  PyPI would have been standing in for. If a PyPI presence is wanted
  later, note that `corpus` (abandoned since April 2018) and
  `corpus-mcp` are both taken; `corpuscle`, `taxon-corpus`,
  `biodiversity-corpus`, `corpus-literature` and similar are free, and a
  distribution rename costs nothing user-visible — the `corpus` command
  and the `pipeline`/`mcpsrv`/`bib` imports are unaffected.
  **Stated** in INSTALL.md §"How corpus is distributed", framed as a
  decision rather than an omission, and the release-pinned install line
  is bumped to `@v1.0.0` — plus a new CONTRIBUTING.md release step to
  bump it, since it names a tag that does not exist until the tag step
  and had therefore already drifted twice.
- [x] **Pull #177 into 1.0** (decided at execution, 2026-08-05). BibTeX
  escapes were served verbatim — five titles in a 702-paper build read
  `A.Braun \& Vatke`. Decoding happens at the `bib/parser.py` boundary
  (escaped specials, group-protection braces, accent commands), with the
  inverse escaping added to `corpus bib export` because a bare `%` in an
  exported `.bib` comments out the rest of the line. The remaining
  post-review issues (#164–#176, #178–#180) defer to v1.1.
- [ ] **Release ritual** per CONTRIBUTING.md §Releasing, including
  **step 8** (prune the roadmap, open the next section) — the step
  skipped at v0.6.0, which is why this document described a completed
  cycle in the present tense for two months.
  **Step 6 is deliberately skipped:** `origin/v0` already exists but
  points at `v0.3.0`+1 (`57fd914`), and v0.6.0 is tagged and immutable,
  so a maintenance branch buys nothing until a v0 backport is actually
  planned. Revisit #48 if one is.

### Carry (tracked, not committed this cycle)

- [ ] **Move the docling pin forward**
  ([#98](https://github.com/caseywdunn/corpus/issues/98) follow-up).
  Still `docling==2.94.0`. Reproduce on an arm64 Mac, determine whether
  2.95/2.96 broke MPS extraction via an API change or an upstream bug,
  then advance deliberately. Needs Apple-Silicon hardware.
- [ ] **Reference reconciliation**
  ([#155](https://github.com/caseywdunn/corpus/issues/155)).
  `get_missing_references` is dominated by unreconciled citation-string
  variants of works already in the corpus —
  `resolve_reference "Bigelow 1911 The Siphonophorae"` returns 53
  matches, one canonical and ~50 variants. Fixing it properly means DOI
  normalization, block-and-cluster canonicalization, alternate-DOI
  aliasing, junk filtering, and probably an auditable LLM adjudication
  pass: a cycle of its own, not a 1.0 item. **Consider the tool-level
  mitigation (§6) as a 1.0 sliver** — collapse candidates by normalized
  (author, year, title) before ranking, so the same work stops appearing
  as N rows. That makes the tool honest without fixing reconciliation,
  and it is additive.

### Deferred to v1.1 by decision (2026-08-05)

Filed after the waves closed — from the pre-release UX review
(#164–#173) and the viburnum build (#174–#177) — and triaged at the
release rather than absorbed into it. Only **#177** was pulled in; the
rest wait. None blocks an install, which is the bar this cycle set.

The freeze permits all of them: every one is a bug fix or an additive
change under [API_STABILITY.md](API_STABILITY.md), so none needs to wait
for a major.

- **Correctness on the served surface.** #175 (taxonomic-authority
  linking assumes zoological authorship, so `get_original_description`
  is structurally dead for any botanical corpus — 889/913 viburnum taxa
  had authorship, 0 with a year), #165 (lexicon translations match only
  uninflected forms, zeroing anatomy coverage on German papers —
  Eschscholtz prints `Luftblasen`, the lexicon has `Luftblase`), #166
  (a hub work's depth-1 `get_citation_graph` payload can exceed MCP
  transport limits while `truncated: false` stays accurate), #155
  (reference reconciliation, above).
- **Extraction quality.** #176 (Tesseract OSD overrides a p=1.0
  language detection with no union or `eng` fallback, and there is no
  per-document override — 6 of 123 OCR'd papers), #164 (abbreviated
  genus binomials, `Ph. pelagica`), #172 (`visual_script` is never
  populated, leaving the text-layer/visual cross-check inert).
- **Operator surface.** #170 (progress heartbeat during long per-document
  stages), #169 (`--filter-gate` silently ignored without
  `--list-hashes` — the hint was fixed, the underlying flag was not),
  #168 (surface the naive-chunker fallback in `corpus status`), #174
  (extend per-stage fingerprints to bib entries, filenames and config
  keys, so editing the input bib stops being a silent no-op).
- **Housekeeping.** #173 (`pipeline/CITATION.cff` is a hand-maintained
  copy with nothing enforcing sync — currently byte-identical, so this is
  a guard rather than a fix), #171 (README ambiguous about where
  `instructions.md` lives).
- **Features, additive by construction.** #178–#180 (a `skills/` plugin
  directory, a clade-monograph skill, and a README quick start).
- **[PR #144](https://github.com/caseywdunn/corpus/pull/144)** from
  @beroe — original-description linking against in-corpus works
  (423 → 505 of 598 ctenophore taxa). Mergeable and substantive, but it
  is a 4-commit external change that has never run CI, and it touches
  `bib/authority.py`, which #154 rewrote during this cycle. It wants a
  review and a CI run, not a last-step merge. Note it overlaps #175:
  both are `parse_authority`, from opposite ends.

### Scope flag (decide at execution)

**#159 is the largest item and the least like the rest.** Waves 0–2 are
pins, one-liners, docs, and bug fixes; `corpus prefetch` is a new
command with a new CLI surface. It earns its place — model downloads are
the single most likely first-run failure for a user on a restricted
network, and #153 and #139 are both that failure in different clothes —
but if the cycle needs to shrink, the trim is to ship the `corpus check`
probe and the `HF_HOME` documentation without the command. The probe is
what converts a mid-run failure into a pre-flight one; the command is
convenience on top.

**#143 is in scope by decision.** It is a query-surface correctness gap
on the surface 1.0 freezes, which argues for now rather than later, but
no external report has it biting anyone. It is the second trim point.

### Landmines

- **#156 blocks verification of everything else.** No fix in this cycle
  can be validated on a clean install until the `mcp` pin lands. Do not
  parallelize Wave 1 work against an env that carries a hand-installed
  `mcp==1.29.0` and assume the result generalizes.
- **Pinning all nine deps (#158b) can fail the solve.** Pinning
  `environment.yaml` and `requirements.txt` to today's resolved versions
  may surface a conflict the cached env was papering over. Do it on a
  branch with T3 available, not in the release commit.
- **#154 §1 and §3 must land together.** §3 tightens strict
  enforcement; §1 stops leaking advisory metadata into permissive use.
  Landing §1 alone widens the `get_figure_roi_image` hole in relative
  terms; landing §3 alone leaves ordinary `report` use still
  self-censoring on 86% of the corpus.
- **#154 §1 narrows a response shape.** Suppressing
  `publishable`/`license_source` under `report` removes fields a v0.6
  client may read. It is the right shape and this is the last cheap
  moment, but it needs a CHANGELOG migration note — and it is the reason
  #154 cannot slip past 1.0.
- **The v0.6 profile landmine still applies, and was only half
  honored.** The v0.6 plan warned: *"The `get_figure_url` → HTTP-fetch
  path must carry the resolved profile … or a strict client leaks
  figures through the unprofiled URL."* `figure_http.py` does check the
  profile — but `get_figure_roi_image` was never considered, which is
  precisely #154 §3. Carry the landmine forward: **every** path that
  returns figure pixels is an enforcement point, and there are now four
  (`figures.py:727`, `:838`, `:594`, `figure_http.py:119`).
- **#142 changes stored provenance.** Stamping unchanged entries flips
  works out of the warning tier, which changes `format_citations` output
  for existing corpuscles on the next import. Desirable, but it is a
  visible behavior change on a frozen tool — changelog it explicitly.
- **#140 re-runs the figure stage** on every existing corpuscle (the
  per-stage resume sees a new stage record). One-time cost, but call it
  out.

### Verification

- **Baselines on `dev` @ `1dbb69a`** — compare against these to spot a
  regression:
  - T0 (`pytest -m "not corpus_required and not resume_scenario"`):
    **623 passed, 2 skipped**
  - `pytest -m corpus_required`: **163 passed, 14 skipped, 4 xfailed** —
    but only *with* the three `--deselect` flags T1 passes
    (`integration.yml:235-237`). Plain `pytest -m corpus_required`
    yields 5 failures and that is expected, not a regression: those
    three functions compare bib-canonical metadata against PDF text and
    legitimately fail on Stepanjants 1970 (Latin title vs Cyrillic OCR),
    Dunn 2005 (`claudanielis ,` spacing), and papers citing no other
    demo paper.
  - SSE smoke (`tools/smoke_test_sse.py demo/output/_serve`): 3 layers,
    **38 tools**
  - Round-2 resume (#71): `embedded=1, skipped=4, failed=0`
- **Unit (T0):** #141 needs a fixture with an escaped brace mid-file
  asserting *all* entries parse and that the skipped-entry count is
  reported; #142 needs an unchanged-`.bib` re-import asserting the
  `bib` tier survives, plus a post-reconcile assertion for #152; #154
  needs a `get_figure_roi_image` refusal under `manuscript` (the §3
  regression), absence of `publishable` under `report` (§1), and a
  five-state `derive_publishable` table (§2); #143 needs
  synonym→canonical, canonical→synonym, and synonym→sibling-synonym
  retrieval; #163 is a metadata assertion. Extend
  `tests/test_freeze_contract.py` to enumerate **every** figure-pixel
  path and assert each one consults the gate — the check that would have
  caught #154 §3.
- **Integration (T1/T2):** the new compose job (#161) asserts
  `docker inspect` health, which is the standing guard for #157. #160's
  pre-flight should be asserted in the negative too — a deliberately
  emptied tessdata dir must make `corpus check` warn.
- **Clean-room (T3):** the §0 gate —
  [`clean-room.yml`](../.github/workflows/clean-room.yml), weekly plus
  `workflow_dispatch`. Dispatch it and require green immediately before
  the tag. It is the only automated tier that exercises a cold model
  download and the real compose file; **T3-bare**
  (`dev_docs/ec2_smoke.sh`) remains the manual pre-release check for the
  bare-host bootstrap.
- Consider adding `tools/` to the pyflakes gate in
  `tests/test_no_undefined_names.py` — it lints `pipeline/`, `mcpsrv/`,
  and `bib/` only, so the four Python scripts under `tools/` never get
  the NameError check #75 built it for.

## 2. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v1.0 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; citation-trust gap closed by [#79](https://github.com/caseywdunn/corpus/issues/79) in v0.5, but provenance preservation is still incomplete — see [#142](https://github.com/caseywdunn/corpus/issues/142) / [#152](https://github.com/caseywdunn/corpus/issues/152) in Wave 2; synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14)) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Indices in place; the corpus-scale vision run + figure-coverage audit is still undone (untracked — see §4) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place; cache cost addressed by dossier tools [#76](https://github.com/caseywdunn/corpus/issues/76) in v0.5 |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Indices in place; figure retrieval is synonym-blind until [#143](https://github.com/caseywdunn/corpus/issues/143) (Wave 3); corpus-scale vision run still undone (untracked — see §4) |

## 3. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual. Note **step 7** — pruning this document at release time
— which was skipped at v0.6.0.

## 4. Carryover to v1.1+ (deferred out of v1.0 by scope)

Not in the installability cycle — either net-new features (safe to add
after 1.0 without breaking the frozen surface) or work awaiting
motivation. Carried so they stay visible.

- **Container distribution image.** The only channel that can ship the
  full native toolchain in one artifact: Docker is already a
  prerequisite for Grobid, `docker-compose.yml` could bring up grobid +
  corpus together, and #153's HPC user already runs Apptainer, which
  pulls straight from a Docker registry. Costs: bind-mounting the PDF
  directory, GPU passthrough for the local VLM, and a large image with
  torch in it. Worth its own issue once 1.0's install path is verified.
- **`export_to_disk` + `suggest_command`**
  ([#88](https://github.com/caseywdunn/corpus/issues/88) Part 2) and the
  **`corpus export` CLI** ([#93](https://github.com/caseywdunn/corpus/issues/93))
  — bulk-export surface for "the LLM isn't the consumer" workflows.
  Additive, so they don't constrain the freeze.
- **`verify_claim`** ([#123](https://github.com/caseywdunn/corpus/issues/123)).
  Per-claim ledger anchoring as a thin similarity-only wrapper over
  `get_chunks_for_topic`. New tool — post-freeze by construction.
- **Drift detection** ([#80](https://github.com/caseywdunn/corpus/issues/80)).
  Pre-run diff of resolved config + per-input content SHAs.
- **Column-store shape for `lexicon_matrix`**
  ([#83](https://github.com/caseywdunn/corpus/issues/83)). Token saving
  at large-matrix scale; held pending a prompt-suite analysis showing it
  matters.
- **Figure-number extraction: the non-caption cases.**
  [#16](https://github.com/caseywdunn/corpus/issues/16) is **closed** —
  it landed the parsing half (`Taf. III.`, `Tab. XII.`, `Plate IV.`,
  Roman→Arabic normalization, 49 tests) of the ~538-of-1,787-paper gap.
  What remains is papers with no caption at all, or a caption not near
  its image, which needs vision OCR or a body-text-mention fallback.
  **Untracked** — file an issue if picked up.
- **Vision pass corpus-scale validation.**
  [#11](https://github.com/caseywdunn/corpus/issues/11) is **closed** as
  carried-out-in-code: `corpus run` invokes the vision pass whenever
  `figures.panel_detection` selects a vision backend and the host
  capability check passes. What remains is the *operational* task — a
  corpus-scale run on Bouchet plus a figure-coverage audit (count
  figures with `pass3c_status` set; sum `missing_figures[]` lengths).
  **Untracked** — release-validation work, not coding work.
- **Evaluate Cloud Run vs the EC2+ALB stack**
  ([#89](https://github.com/caseywdunn/corpus/issues/89)). Deployment
  decision, not part of the MCP API contract.

## 5. Out of scope (longer horizon)

- [#124](https://github.com/caseywdunn/corpus/issues/124) — whether to
  expand the server-side LLM-call surface at all. `translate_chunk` was
  removed in v0.6, leaving the MCP server a purely deterministic
  retrieval layer. Re-opening that is a direction question, not a task.
- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to warrant
  its own plan section when picked up.
- [#13](https://github.com/caseywdunn/corpus/issues/13) — Geographic
  mention layer. Deferred to v2.0+; the mention-layer surface is likely
  to be reworked at the major-version boundary.
- [#38](https://github.com/caseywdunn/corpus/issues/38) — Embedding
  model migration path. Design-only; implement when a model swap is
  actually needed.
- [#39](https://github.com/caseywdunn/corpus/issues/39) — MCP server
  lazy index loading. Premature at 1.8K papers; documented as a known
  scaling cliff; revisit when a corpuscle pushes 10K+.
- [#5](https://github.com/caseywdunn/corpus/issues/5) — Streamable HTTP
  transport with OAuth. Deferred indefinitely. SSE + bearer-token works
  for the ~20-collaborator deploy target.
- Multi-region failover, autoscaling, or blue/green deploys for the AWS
  served bundle. Single-instance per corpuscle is fine until it isn't.
- Authentication beyond bearer tokens / OAuth (Cognito, institutional
  SSO).
- Mirror to Cloudflare R2 / Backblaze B2 for cost. Defer until S3 egress
  shows up on a bill.
- A thin HTML/web UI on top of the MCP server. Out of scope until the
  MCP-only experience has actual non-Claude-Desktop users.

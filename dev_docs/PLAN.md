# PLAN.md — Corpus pipeline (v0.5)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing
→ MCP-serving stack. v0.2 (2026-05-08) hardened internals: per-stage
resume, input fingerprints, structured failure schema, quality gates,
and an ALB-fronted multi-organism deploy. v0.3 (2026-05-11) collapsed
the user surface: one CLI installed once, a per-corpuscle
`config.yaml`, automatic database management, capability-aware
defaults. v0.4 (2026-05-17) closed the silent-failure modes external
testers and clean-room HPC hosts surfaced, and gated the next cycle
on tiered CI that didn't exist before.

**v0.5 is the served-bundle quality cycle.** v0.4 made the pipeline
that *produces* a corpus reliable. v0.5 turns to the surface a
downstream LLM client *consumes* — the MCP tool surface — and closes
the two trust-breaking gaps external evaluators surfaced there:
amalgamated citations the model assembles in its own context window,
and an MCP access pattern whose `for each X, call get_Y` shape burns
~2 M cache_read tokens on enumerative prompts. Both came from real
evaluator feedback against the v0.3 / v0.4 build.

The cycle is non-breaking on the pipeline side — v0.4.x corpuscles
build as-is — but adds MCP tools and reshapes how clients are
expected to consume the bibliography surface. Operators upgrading
their served bundle pick up the new surface; client-side prompts that
hand-wrote citations evolve to route through `format_citation`.

This is *not* v1.0. 1.0 carries an API-stability signal, and the
citation-trust gap that motivates v0.5 has to be solved and proven
in production before that commitment is on the table. Save 1.0 for
after v0.5 runs a cycle.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](DEPLOY.md); the
full v0.4 plan archived in the
[v0.4.0 tag's history](https://github.com/caseywdunn/corpus/blob/v0.4.0/dev_docs/PLAN.md).
Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## 1. v0.5 punch list

### Citation grounding — stop LLM-side hallucination

External evaluators have caught hallucinated references in
documents generated against the v0.3 / v0.4 served bundle —
amalgamated entries with the right title but the wrong author
and journal. One incident: a taxonomist found his own paper in a
generated reference list attributed to a different author and a
different journal. Citations are the poster-child for "model
misbehaving"; a single such miss disproportionately shreds trust.

Root cause is not the corpus data. `biblio_authority.sqlite` is
authoritative — the reconciliation cascade plus user-edited `.bib`
import put a clean, provenance-tagged record in front of the model
(reconciliation defects #1 and #2 are closed). Hallucination happens
client-side, because today's MCP surface returns structured fields
and leaves the LLM to recombine them into citation prose in its own
context window. That recombination is where amalgamation happens.

- [x] **`format_citation` MCP tool + provenance cascade**
  ([#79](https://github.com/caseywdunn/corpus/issues/79)).
  Shipped in four commits on dev: schema + importer stamp
  (219a4a2), `bib/format.py` with the author-year formatter
  (4d16064), `BiblioAuthority.provenance()` tier classifier
  (fc3a7e9), and the `format_citation` MCP tool itself wiring
  the three together (e8f075a). Public warning footnotes per
  tier match the issue spec verbatim. `vancouver` / `bibtex`
  styles + the `bibliography.citation_style` config field are
  deferred to follow-ups — only one style currently exists in
  `SUPPORTED_STYLES`, so a config field with no observable
  effect would be half-functional cruft. Will wire when a second
  style lands.
- [x] **`mcpsrv/default_instructions.md` rewrite**
  (part of [#79](https://github.com/caseywdunn/corpus/issues/79)).
  "Defer to the corpus bibliography" section rewritten as an
  explicit routing rule: call `format_citation`, paste its
  `formatted` and `inline` verbatim, preserve `warning` verbatim,
  and handle `not_found` / `ambiguous` error shapes without
  fabricating. Lands in every session's context. Implementation
  in the next commit after PLAN tick.
- [x] **Citation-grounding quality tests**
  (part of [#79](https://github.com/caseywdunn/corpus/issues/79)).
  Harness shipped in `tests/test_prompt_quality.py`: real Claude
  API round-trips with `format_citation` wired as a callable tool,
  scoring on expected-citation emission + parenthetical-token
  trace-back + forbidden-hallucination co-occurrence windows.
  Gated behind `RUN_PROMPT_QUALITY=1` + `ANTHROPIC_API_KEY` so
  default CI doesn't burn API budget. Schema for the new
  `prompts:` block in ground-truth YAMLs documented in
  `dev_docs/TESTING.md`. Scoring functions have in-process unit
  tests that always run. First concrete prompt panels (and CI
  integration that sets the API key secret in the release-time
  lane) are tracked as follow-ups — needs the demo to be built
  with a current biblio_authority.sqlite to look up real work_ids.

### Cache-friendly MCP dossier tools

External-evaluation feedback on token / cache cost: the current
`for each X, call get_Y` patterns cost ~1.95 M cache_read tokens
on enumerative prompts (e.g. p02 adiastola anatomy), because each
`get_chunks_for_taxon` lands ~5–10 k tokens of full text that every
subsequent turn re-reads. Pre-aggregated dossier tools return
structured indices + IDs + headings, with full text pulled
selectively. Worked example on p02 projects ~85 % token / ~92 %
cache_read reduction. Pairs naturally with #79: both reshape the
MCP surface for the way LLM clients actually consume it.

- [ ] **Dossier tool suite**
  ([#76](https://github.com/caseywdunn/corpus/issues/76)).
  Six new tools per the issue's priority order:
  `corpus_summary`, `get_taxon_dossier`, the
  `get_chunks(ids)` drill-down pair, `get_taxon_lexicon_slice`,
  `lexicon_matrix`, `get_figure_dossier_*`, and batched
  `get_papers`. 7–10 days across the six; the underlying indexes
  are all on disk already, so v0.4's hardening doesn't block it.
  Design lives in #76.

### Quick fixes + carryover

Smaller items that fit alongside the two big-ticket pieces.

- [x] **`taxonomy.root_id` rejects DwC-A LSID format**
  ([#78](https://github.com/caseywdunn/corpus/issues/78)).
  Schema relaxed to `Optional[Union[int, str]]`; worms keeps the
  integer AphiaID path, dwc/dwca now accept either a bare integer
  taxonID or an LSID string. Downstream code already wrapped with
  `str(...)`, so no callers changed. Regression tests pin both
  shapes. Landed on dev (7be13f8); closes on release merge.
- [x] **Install docs / conda prescription**
  ([#77](https://github.com/caseywdunn/corpus/issues/77)).
  README + INSTALL.md reframed: the requirement is an
  arm64-native conda (Anaconda, Miniconda, or Miniforge all
  qualify), not miniforge specifically. New `conda info | grep
  platform` pre-gate; `corpus check` now hard-fails on Rosetta'd
  Python on macOS (post-gate). `CONDA_SUBDIR=osx-arm64` is
  documented as the alternative for users who don't want a second
  conda distribution. Installation section restructured so
  Supported platforms comes first, macOS-specific items are
  consolidated in one subsection, and Grobid startup moved out to
  its own h2 (it's a run-time, not install-time, concern).
  Landed on dev (b92e355 + follow-ups); closes on release merge.
- [ ] **Drift detection.** Hash the resolved config (input paths +
  per-input content SHA) into a corpuscle-side state file so a
  re-run that's doing more than expected can show *why*. Carried
  from v0.3 → v0.4; not yet filed.
- [x] **Success summary on clean `corpus run`**
  ([#57](https://github.com/caseywdunn/corpus/issues/57)).
  New `_write_run_log` helper in `pipeline/cli.py` writes the
  same content as `corpus status --report` to
  `<output_dir>/run.log` at the tail of `corpus run`. Trailing
  success message points at the file. Overwrites prior log; dry-run
  short-circuits before writing. Landed on dev (d365792); closes
  on release merge.
- [x] **Lexicon YAML loader bug.** Tightened: `load_lexicon`
  now rejects list-shaped term entries with a ValueError that
  names the category, term, offending type, and corrected shape;
  `main.py` narrowed its swallow to file / yaml errors so schema
  violations propagate instead of silently degrading to no-op
  annotation. Landed on dev (f4b36b0).

## 2. Carryover to v0.6

Out of v0.5 by scope, carried so they stay visible.

- **Figure-number extraction on old/scanned papers**
  ([#16](https://github.com/caseywdunn/corpus/issues/16)).
  ~538 of 1,787 papers have unparsed figure numbers; modest
  heuristic work.
- **Vision pass corpus-scale validation**
  ([#11](https://github.com/caseywdunn/corpus/issues/11)).
  Operational run + figure-coverage audit on Bouchet against the
  v0.5 build. Release-validation, not coding.

## 3. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v0.5 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; citation-trust gap addressed by [#79](https://github.com/caseywdunn/corpus/issues/79); synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14), v0.6+ candidate) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place; cache cost addressed by dossier tools [#76](https://github.com/caseywdunn/corpus/issues/76) |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |

## 4. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual.

**v0.5 release note:** non-breaking on the pipeline side; v0.4.x
corpuscles build as-is. The MCP surface gains tools (`format_citation`,
six dossier tools) but does not remove or rename any existing tools —
clients that ignore the new surface continue to work. Operators who
*want* the citation-grounding behavior reshape their prompts (or rely
on the updated `default_instructions.md`) to route through
`format_citation`.

## 5. Out of scope for v0.5

- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to
  warrant its own plan section when picked up; v0.6+ candidate.
- [#13](https://github.com/caseywdunn/corpus/issues/13) — Geographic
  mention layer. Deferred to v2.0+; the mention-layer surface is
  likely to be reworked at the major-version boundary.
- [#38](https://github.com/caseywdunn/corpus/issues/38) — Embedding
  model migration path. Design-only; implement when a model swap is
  actually needed (BGE-M3 v2 or a switch to a different family).
- [#39](https://github.com/caseywdunn/corpus/issues/39) — MCP server
  lazy index loading. Premature at 1.8K papers; documented as a
  known scaling cliff; revisit when a corpuscle pushes 10K+.
- [#5](https://github.com/caseywdunn/corpus/issues/5) — Streamable
  HTTP transport with OAuth. Deferred indefinitely. SSE +
  bearer-token works for the ~20-collaborator deploy target.
- Multi-region failover, autoscaling, or blue/green deploys for the
  AWS served bundle. Single-instance per corpuscle is fine until
  it isn't.
- Authentication beyond bearer tokens / OAuth (Cognito,
  institutional SSO).
- Mirror to Cloudflare R2 / Backblaze B2 for cost. Defer until S3
  egress shows up on a bill.
- A thin HTML/web UI on top of the MCP server. Out of scope until
  the MCP-only experience has actual non-Claude-Desktop users.

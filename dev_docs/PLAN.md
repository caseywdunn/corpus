# PLAN.md — Corpus pipeline (v0.6)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing
→ MCP-serving stack. v0.2 (2026-05-08) hardened internals: per-stage
resume, input fingerprints, structured failure schema, quality gates,
and an ALB-fronted multi-organism deploy. v0.3 (2026-05-11) collapsed
the user surface: one CLI installed once, a per-corpuscle
`config.yaml`, automatic database management, capability-aware
defaults. v0.4 (2026-05-17) closed the silent-failure modes external
testers and clean-room HPC hosts surfaced, and gated the next cycle
on tiered CI. v0.5 (2026-05-29) was the served-bundle quality cycle:
it closed the two trust-breaking gaps external evaluators found in the
MCP surface — LLM-side citation amalgamation (`format_citation` +
provenance cascade, [#79](https://github.com/caseywdunn/corpus/issues/79))
and the ~2 M cache_read enumerative access pattern (the six-tool
dossier suite, [#76](https://github.com/caseywdunn/corpus/issues/76)) —
plus token-efficiency follow-ups
([#81](https://github.com/caseywdunn/corpus/issues/81),
[#82](https://github.com/caseywdunn/corpus/issues/82),
[#84](https://github.com/caseywdunn/corpus/issues/84),
[#86](https://github.com/caseywdunn/corpus/issues/86)). The release
itself surfaced one more silent-failure mode and its cause: an
unpinned docling/torch stack let a same-day docling release silently
break macOS arm64 extraction, which shipped an empty bundle reporting
success. v0.5 closed both — pinning the ML stack to a known-good set
([#98](https://github.com/caseywdunn/corpus/issues/98)) and making
`corpus run` fail loudly on a zero-yield extraction
([#99](https://github.com/caseywdunn/corpus/issues/99)).

**v0.6 is the road-to-1.0 cycle.** 1.0 carries an API-stability
signal, so the through-line is *no new feature tools* — instead,
finalize the public MCP tool surface, fix known correctness bugs, and
harden the production/operations experience. The organizing principle:
**any change to an MCP tool's default, signature, or response shape
must land in v0.6**, because after 1.0 those become breaking changes.
The spine of the cycle is therefore a one-time "API-freeze pass" over
the tool surface (Group A), with a focused set of bug fixes (Group B)
and ops polish (Group C) alongside. Writing the formal 1.0
API-stability *policy* doc is deliberately left for the 1.0 release
itself; v0.6 is code and fixes.

> **Planning note — #88 Part 1 is not merged.** The
> [#88](https://github.com/caseywdunn/corpus/issues/88) body describes
> "Part 1" changes as implemented "in this branch," but `git log --all`
> shows no trace and the tools confirm it: `format_citations` (plural)
> does not exist, `with_cites`/`parent_chain` are absent, and
> `lexicon_matrix` still returns the full grid. So the Group-A Part-1
> items below are *net-new implementation*, not settling existing code.
> See the **scope flag** under §1.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](DEPLOY.md); the
full v0.5 plan archived in the
[v0.5.0 tag's history](https://github.com/caseywdunn/corpus/blob/v0.5.0/dev_docs/PLAN.md).
Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## 1. v0.6 punch list

### Group A — MCP tool-surface API-freeze pass

The spine of the cycle. **Sequence matters**: implement all signature
changes first, remove redundant tools second, reconcile naming /
error-shape last (the reconciliation is the freeze gate — doing it
earlier guarantees a second rename pass). Phase 1 items are independent
and parallelizable across issue branches; Phase 2 depends on Phase 1;
Phase 3 is last.

**Phase 0 — instrumentation scaffold** (do first; unblocks #91/#90 and
the uniform error shape)

- [ ] **Per-tool instrumentation shim** in
  [mcpsrv/app.py](../mcpsrv/app.py) at registration (tool name, args
  summary, latency, ok/error → module counter + logger). Internal,
  additive. Feeds #91 structured logging, #90 success/failure counts,
  and the Phase-3 uniform error-row shape. Reuse the existing
  `_validated_limit()` (`app.py:53`) as the `limit` standardization
  lever.

**Phase 1 — per-tool signature changes**

- [ ] **`lexicon_matrix` summary-by-default**
  (part of [#88](https://github.com/caseywdunn/corpus/issues/88) Part 1;
  `mcpsrv/tools/lexicon.py:89`). Add `detail: bool=False`; default
  returns per-term totals, opt-in returns the full grid. **Breaking
  default** — fixes the 71 MB → 684 KB runaway.
- [ ] **Batched `format_citations`**
  ([#88](https://github.com/caseywdunn/corpus/issues/88);
  `mcpsrv/tools/bibliography.py`). New `format_citations(queries |
  work_ids | paper_hashes, style)` reusing the existing single-resolve
  body. Must exist before the singular is removed in Phase 2.
- [ ] **`with_cites=True` on `get_chunks_for_topic`**
  ([#88](https://github.com/caseywdunn/corpus/issues/88);
  `mcpsrv/tools/chunks.py:243`). Attach in-text cite refs per parent
  paper. **Breaking default.** *(Enhancement-leaning — trimmable; see
  scope flag.)*
- [ ] **`search_taxon` `parent_chain`**
  ([#88](https://github.com/caseywdunn/corpus/issues/88);
  `mcpsrv/tools/taxonomy.py:40`). Additive ancestry walk.
  *(Enhancement-leaning — trimmable; see scope flag.)*
- [ ] **Caption preview default on older figure tools**
  ([#85](https://github.com/caseywdunn/corpus/issues/85);
  `mcpsrv/tools/figures.py`). In `get_figures_for_taxon` /
  `get_figures_for_lexicon_term` (full caption at `:156`, `:240`) swap
  `caption_text` for `_caption_preview(...)` (helper already at `:257`)
  and add `full_caption: bool=False`. **Breaking** (default field
  truncated).
- [ ] **Breadth + edge caps on `get_citation_graph`**
  ([#87](https://github.com/caseywdunn/corpus/issues/87);
  `mcpsrv/tools/bibliography.py:205-276`). Add `max_edges_per_node` /
  `max_total_edges`, rank each node by `cited_by_count` in
  `_walk_citations` (`:258`), add a `truncated` field. Keep defaults
  generous → additive.
- [ ] **Output-type profile ontology — per-call `profile=` selection**
  ([#101](https://github.com/caseywdunn/corpus/issues/101), supersedes
  [#94](https://github.com/caseywdunn/corpus/issues/94)). Land **on the
  same branch as #85** (both touch `figures.py`). Output type is a
  **client/session property carried per call, not a corpus/server
  attribute** — the shared SSE deploy (~20 clients, one index) means a
  server-start flag or per-corpuscle-config default can't distinguish a
  user's concurrent chat / internal-report / manuscript sessions. Add an
  optional `profile=` arg (global built-in vocabulary `report` /
  `manuscript` / `presentation`) to the gated figure + citation tools;
  remove `--allow-unpublishable` (`mcpsrv/main.py:125-130/252`) for a
  `--default-profile` server *fallback* (default `report`); enforce
  `figure_licensing` at `figures.py:680,786` + `figure_http.py:103`; add
  `list_output_profiles` / `get_active_profile` discovery tools.
  **Breaking** (CLI flag removed; permissive fallback; figure response
  shape gains `license`/`attribution` fields). Startup warning when the
  fallback is permissive — security-relevant default. **Freeze-critical
  core** = the `profile=` param + stable profile vocabulary + figure gate
  (these change signatures/defaults, so must land pre-1.0); the further
  axes (`require_attribution` emission, `citation_provenance: strict`,
  `excerpt_max_words`) are additive and may follow post-1.0 —
  `citation_provenance: strict` depends on #100.

**Phase 2 — tool removals** (after Phase 1)

- [ ] **Remove the three redundant singular tools**
  ([#88](https://github.com/caseywdunn/corpus/issues/88) §2.3). Delete
  `format_citation` (`bibliography.py:352`), `get_paper`
  (`papers.py:208`), `get_chunk` (`papers.py:297`) — the plurals
  `get_papers` / `get_chunks` already cover the singleton case.
  **Breaking (tools removed).** Rewire
  [mcpsrv/default_instructions.md](../mcpsrv/default_instructions.md)
  (the citation routing rule) and `dev_docs/MCP_TOOLS.md`. Do **not**
  touch `bib.format.format_citation` — that is a different symbol (the
  pure string formatter); leave its import intact.

**Phase 3 — consistency reconciliation** (last; the freeze gate)

- [ ] **Reconcile naming + error shape across the tool surface.** Pick
  one convention — `limit` everywhere (rename figure-tool
  `max_figures` / `max_linked_chunks`, or document the dossier
  exception) — and apply the Phase-0 uniform error-row shape
  (`{"error", "code"}`) to every `return [{"error": ...}]` site. One
  consolidated CHANGELOG entry + a migration table. This is the surface
  that 1.0 freezes.

### Group B — correctness bugs

Runs in parallel with Group A (different subtree).

- [ ] **Bib-provenance preservation through import + reconciliation**
  ([#100](https://github.com/caseywdunn/corpus/issues/100), follow-up to
  [#79](https://github.com/caseywdunn/corpus/issues/79)). Over-zealous
  `format_citation` warnings (reported by shchurch): bib-imported
  references all warn because the `bib` provenance tier isn't stamped or
  preserved. Three causes — (a) `bib_imported_at` is only stamped on
  *changed* rows (`bib/importer.py:230-232`), so an unchanged re-import of
  the authoritative `.bib` leaves every row warned; (b)
  `_merge_phase1_into_ghost` (`bib/reconcile.py:313-325`) keeps the
  GROBID-derived ghost row and deletes the bib-stamped phase-1 row without
  carrying `bib_imported_at`/bib fields forward, so reconciliation against
  corpus-paper bibliographies drops bib authority; (c)
  `find_matching_work_id` (`bib/importer.py:80-111`) matches only on
  `corpus_hash`/`doi`. **Invariant**: a reference present in the
  user-edited `.bib` stays authoritative `bib` provenance, including after
  reconciliation. Unblocks #101's `citation_provenance: strict` axis.
- [ ] **`build_taxon_mentions` freshness**
  ([#95](https://github.com/caseywdunn/corpus/issues/95)). The mtime
  gate (`pipeline/taxon_mentions.py:191-201`) is correct in isolation;
  the bug is cross-component — a taxonomy refresh can leave per-paper
  `taxa.json` mtime stale while the backbone changed. **Make staleness
  fingerprint-aware**: compare each paper's `pipeline_state.json`
  taxa-stage `input_fingerprint` against the current `taxonomy.sqlite`
  sha256 (`pipeline/main.py:338`), not just mtime. Validate against the
  exact repro (rebuild taxonomy.sqlite → re-run → assert mentions
  re-ingest).
- [ ] **HuggingFace token warning**
  ([#97](https://github.com/caseywdunn/corpus/issues/97)). One shared
  `_hf_env()` setter (`HF_HUB_DISABLE_IMPLICIT_TOKEN=1`) imported by the
  three model-load sites (`pipeline/embeddings.py:183`,
  `pipeline/vision.py:429/434`, `pipeline/extract.py:59`), or set once
  at pipeline entry. Not a contract change.

### Group C — production hardening / UX polish

- [ ] **MCP health checks + refuse-to-serve on degraded capability**
  ([#91](https://github.com/caseywdunn/corpus/issues/91)). Turn the
  silent `(None, None)` degradation in `get_topic_searcher`
  (`mcpsrv/indexes.py:91-132`) into a tracked degraded state; make
  `/healthz` (`mcpsrv/transport.py:32-59`) return JSON capability flags
  + non-200 when a required capability is down; raise hard MCP errors
  (not empty rows) when a tool's backing index is degraded. Consume the
  Phase-0 instrumentation for structured request logging.
- [ ] **Central per-invocation run log**
  ([#90](https://github.com/caseywdunn/corpus/issues/90)). Extend
  `_write_run_log` (`pipeline/cli.py:504`) to write
  `<output_dir>/runs/<timestamp>/run.log` (keep the top-level `run.log`
  as the "latest" copy for #57 back-compat) with argv, resolved config,
  tool/dep versions (`pipeline/version.py`), and success/failure counts.
- [ ] **`corpus debug-pdf` single-PDF debug runner**
  ([#92](https://github.com/caseywdunn/corpus/issues/92)). New subparser
  in the `pipeline/cli.py` router (mirror `run_p` / `status_p`) that runs
  `pipeline/runner.py`'s per-paper pipeline on one PDF with verbose stage
  tracing, no bundle. Update fish completions.
- [ ] **Document the WoRMS `isMarine=0` gap**
  ([#96](https://github.com/caseywdunn/corpus/issues/96)). Docs only:
  the WoRMS DwC-A backbone excludes `isMarine=0` records, so
  freshwater/terrestrial taxa hit silent resolution failures. Add to
  README/INSTALL or a `dev_docs/` taxonomy note.

### Carry (tracked, not committed this cycle)

- [ ] **Move the docling pin forward**
  ([#98](https://github.com/caseywdunn/corpus/issues/98) follow-up).
  Keep `docling==2.94.0`; reproduce on an arm64 Mac, determine whether
  docling 2.95/2.96 broke MPS extraction via an API change to adapt to
  or an upstream bug to report, then advance the pin deliberately.
  Needs Apple-Silicon hardware to verify — not blocking the cycle.

### Scope flag (decide at execution)

#88 Part 1 being unimplemented enlarges Group A beyond what "settle
defaults" implied. Two Part-1 items are enhancement-leaning rather than
fix/consolidation: `with_cites` on `get_chunks_for_topic` and
`search_taxon` `parent_chain`. The others are load-bearing for the
freeze (`lexicon_matrix` summary fixes a runaway; `format_citations` is
required before the singular can be removed). Recommendation: keep all
of Part 1 — these are the right 1.0 default shapes and the #88 eval data
backs them — but this is the natural trim point if the cycle needs to
shrink.

### Landmines

- **Tool removal (#88 §2.3) is highest-risk.**
  `default_instructions.md` tells the served LLM to call
  `format_citation` for *every* citation — it must be rewired to
  `format_citations` or the assistant calls a dead tool.
  `tests/test_format_citation_tool.py` and `tests/test_prompt_quality.py`
  (a real-model eval needing `ANTHROPIC_API_KEY`) are pinned to the
  singular and must be rewritten + re-run. No internal Python caller uses
  the singular MCP tools.
- **#101 is a security-relevant default flip** (permissive figure serving
  when a call omits `profile=`). A public deploy relying on today's
  default-deny becomes permissive on upgrade, and `--allow-unpublishable`
  is removed — prominent CHANGELOG migration note + startup warning. The
  `get_figure_url` → HTTP-fetch path must carry the resolved profile (the
  handed-out URL is enforced independently at `figure_http.py:103`), or a
  strict client leaks figures through the unprofiled URL.
- **Never remove `format_citation` before `format_citations` exists**
  (Phase 1 → Phase 2) or the server has zero citation tool.

### Verification

- **Unit (T0):** rewrite `test_format_citation_tool.py` →
  `format_citations` + a regression test asserting the singular tools are
  no longer registered; extend `test_lexicon_tools.py` (summary vs
  `detail=True`), `test_chunks_for_topic_with_text.py` (`with_cites`),
  `test_figure_dossier.py` (caption truncation + `full_caption`); new
  `test_citation_graph_caps`; per-profile figure gate (#101) including two
  concurrent SSE clients passing different `profile=` and getting
  independent results (the core per-session regression) + attribution
  emitted under `manuscript`; bib-provenance preservation (#100: an
  unchanged `.bib` re-import and a post-import reconcile both keep the
  `bib` tier with no warning); fingerprint-repro
  in `test_taxon_mentions_fast_path.py`, HF-warning-absent,
  `test_healthz`, `test_debug_pdf`, `test_run_log.py` (runs/<ts>/). Add a
  Phase-3 "freeze contract" meta-test enumerating registered tools and
  asserting uniform error shape + `limit` naming.
- **End-to-end (T1 Linux / T2 macOS):** demo build + MCP SSE round-trip
  exercises #91 healthz, #90 run.log, #92 debug-pdf, the #101 profile gate
  (per-call `profile=` over a shared SSE server), and the
  §2.3 tool-list. The §2.3 rename **mandates** an `ANTHROPIC_API_KEY` run
  of `test_prompt_quality.py`.
- Re-run the #88 eval suite (`eval/run_suite.py`, lite suite) after the
  surface changes to confirm the turn/volume wins and no accuracy
  regression.

## 2. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v0.6 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; citation-trust gap closed by [#79](https://github.com/caseywdunn/corpus/issues/79) in v0.5 (provenance-preservation follow-up [#100](https://github.com/caseywdunn/corpus/issues/100)); synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14)) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place; cache cost addressed by dossier tools [#76](https://github.com/caseywdunn/corpus/issues/76) in v0.5 |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |

## 3. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual.

## 4. Carryover to v1.1+ (deferred out of v0.6 by scope)

Not in the road-to-1.0 freeze cycle — either net-new features (safe to
add after 1.0 without breaking the frozen surface) or work awaiting
motivation. Carried so they stay visible.

- **`export_to_disk` + `suggest_command`**
  ([#88](https://github.com/caseywdunn/corpus/issues/88) Part 2) and the
  **`corpus export` CLI** ([#93](https://github.com/caseywdunn/corpus/issues/93))
  — bulk-export surface for "the LLM isn't the consumer" workflows. New
  tools; additive, so they don't constrain the 1.0 freeze.
- **Drift detection**
  ([#80](https://github.com/caseywdunn/corpus/issues/80)). Pre-run diff
  of resolved config + per-input content SHAs.
- **Column-store shape for `lexicon_matrix`**
  ([#83](https://github.com/caseywdunn/corpus/issues/83)). Token saving
  at large-matrix scale; held pending a prompt-suite analysis that shows
  it matters.
- **Figure-number extraction on old/scanned papers**
  ([#16](https://github.com/caseywdunn/corpus/issues/16)). ~538 of 1,787
  papers have unparsed figure numbers.
- **Vision pass corpus-scale validation**
  ([#11](https://github.com/caseywdunn/corpus/issues/11)). Operational
  run + figure-coverage audit on Bouchet.
- **Evaluate Cloud Run vs the EC2+ALB stack**
  ([#89](https://github.com/caseywdunn/corpus/issues/89)). Deployment
  decision, not part of the MCP API contract.

## 5. Out of scope (longer horizon)

- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to
  warrant its own plan section when picked up.
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

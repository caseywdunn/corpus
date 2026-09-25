# v1.5 correctness release: session handoff

For the separately authorized September 25 rebuild on Bouchet and newly
available local acceptance assets, see [V1_5_BOUCHET_BUILD.md](V1_5_BOUCHET_BUILD.md).
The dated erenna observations below remain historical.

Recorded **2026-09-24, 09:22 EDT / 13:22 UTC**. This is a dated handoff,
not a live build dashboard. Recheck the tracker and build before acting.
The siphonophore paths below describe this validation run, not requirements
for other corpora or machines.

## Start here on another machine

Fetch and check out **`release/v1.5-correctness`**, then read this note,
[PLAN.md](PLAN.md), and [V1_5_SOURCE_ACCEPTANCE.md](V1_5_SOURCE_ACCEPTANCE.md).
The implementation/evidence commit before this handoff is
`384911ad4a569da0f38258c785a682bde3bad47f`, verified pushed to GitHub.
The handoff commit is its documentation-only successor.

**The full candidate build is still running on `erenna`. A Git checkout does
not bring that process, its inputs, model cache, outputs, or raw logs to a new
machine.** Leave it running and obtain its result from that host before
starting another full build. If the host is unavailable, report that limitation;
do not interpret an inaccessible run as a failed or completed run.

The next substantive step is to inspect the completed candidate's phase audits,
then perform the full-corpus bibliography and frozen retrieval acceptance.
Do not restart already completed gold acceptance merely because the session
has moved.

## Scope, branches, and authorization

- v1.5 is the correctness release. Enhancements are preserved for v1.6 on
  `enhancements/v1.6`; it and `dev` remain at
  `44977fea1657985f075f7e3bafc7833d2ab914cd`.
- The temporary integration branch is `release/v1.5-correctness`, based on
  published v1.4.0 (`b1e34ad`). Follow PLAN's temporary integration arrangement
  rather than merging these fixes into `dev` prematurely. Carry fixes forward
  before enhancements resume; do not reset or force-push either branch.
- [PR #335](https://github.com/caseywdunn/corpus/pull/335) is open and **draft**,
  targeting `main`. No release, tag, deployment, or merge to `main` has been
  performed or authorized.
- User instruction: finish remaining originally scoped issues, then added
  in-scope issues, except where a dependency makes an added issue best first.
  Correctness takes priority over rushing to fit a usage allowance. No numeric
  token budget was specified.
- **#340 is explicitly deferred**. Keep the existing 35-document siphonophore
  gold set as the fixture; do not add the large Hydractinia corpuscle as a new
  fixture. The full reference candidate is a separate acceptance build.
- The main worktree was clean before this note. Two old temporary worktrees
  have five untracked copies: all match files in current or earlier pushed
  release history, with no unique implementation to recover. They were left
  untouched. Completed fixes from agent branches are already integrated.

## Remaining actionable issues

GitHub issue state was checked on September 24. Seven issues remain actionable
in this release; implementation completion is distinct from source acceptance.

| Issue | Scope / plane | Remaining acceptance or work |
| --- | --- | --- |
| [#296](https://github.com/caseywdunn/corpus/issues/296) | Original; build | BibTeX provenance and author precedence fixes/regressions landed. Named audited production merge-history review still needs the September 9 v1.4 authority/bundle. |
| [#305](https://github.com/caseywdunn/corpus/issues/305) | Original; build | Both coordinate backends are fixed. Fresh compatible-GPU inference and visual crop checks remain; the original six-case ROI/control manifest is also unavailable. |
| [#314](https://github.com/caseywdunn/corpus/issues/314) | Original; build | All four named source cases pass normal preparation, Grobid, and materialization. Full-corpus citation overlap and missing-work ranking still need review. Audited v1.4 history is unavailable. |
| [#320](https://github.com/caseywdunn/corpus/issues/320) | Original; build and bounded serve/query | Evaluator, frozen source labels, text/context corrections, and full retained baseline are ready. Complete candidate, compare the same populations, investigate misses, and make only supported further changes. Retrieval improvement is **not yet measured**. |
| [#336](https://github.com/caseywdunn/corpus/issues/336) | Added; build | Next-page caption/plate association fix and regressions landed. Actual Porifera PDF/artifacts `21ff756cca59` remain unavailable for source replay. |
| [#337](https://github.com/caseywdunn/corpus/issues/337) | Added; validation tooling, no plane | General SSE smoke implementation, regressions, and live demo/reference checks pass. A non-reference production bundle is still needed. |
| [#342](https://github.com/caseywdunn/corpus/issues/342) | Added; build | Qwen resized-image coordinates are fixed. Shares #305's GPU/crop acceptance, not a separate model rewrite. |

**#298 is open on GitHub but its v1.5 warning/help component is complete.**
Broader multi-root ingestion is intentionally deferred. #303 was closed after
current source acceptance; added #338, #339, and #341 are also closed.
#297, #330, #333, enhancements, and other open repository issues are not newly
in scope. See PLAN for release-wide gates beyond these seven issues, including
the served corpus audit, remaining clean/incremental matrix, operator/migration
documentation, final checks and release metadata.

## Live full-corpus candidate on erenna

At the timestamp above, the supervisor was alive and the extraction log was
updating. **466 document summaries report success**, with no recorded errors,
stage failures, or failed stage timings in either the top-level summary or its
`processing_summary`. There are quality warnings in 363 of those documents;
successful execution is not a claim of scientific fidelity.

Five `gibberish_after_ocr` alerts have severity `error`: Alvarino_Kimbrell1987,
delle Chiaje1841PlatesVolumes6-7, Aurivillius1898, delle Chiaje1822Plates, and
Alvarino1980b. Review their sources and extraction; these alerts are neither
recorded stage failures nor already adjudicated source defects.

The run was processing document 467, `delle Chiaje1841Volume4.pdf`
(`40664a677b70`), with 145 included pages. These counts describe the frozen
September 22 candidate, not the size of the growing library. Elapsed time was
about 43.5 hours. Taxonomy passed; extraction is ongoing; embedding, post, and
bundling have not started. About 234 GiB was free. No reliable completion-time
estimate has been established.

| Item | Location / identity on erenna |
| --- | --- |
| Candidate root | `/tmp/corpus-v15-full-candidate-20260922` |
| Current run directory | `/tmp/corpus-v15-full-candidate-20260922/runs/20260922T174832.798244Z` |
| Live status and log | `run.json` and `extract.log` under that run directory |
| Build outputs | `/tmp/corpus-v15-full-candidate-20260922/output` |
| Driver and audits | `driver.py`, `audit_candidate.py`, and phase receipts under the candidate root/run directory |
| Supervisor | PID `820851` at observation; verify command and start time before trusting a later PID lookup |
| Pinned code worktree | `/tmp/corpus-v15-full-candidate-code-20260922` at `7b4343f7c466d091568132203db169e9126f14ee` |
| Production implementation | Identical to `ab6353f8d4a8ef59b1be314458d32437913635ec`; later commits contain tooling/evidence/docs |
| Environment | `/home/claude/miniforge3/envs/corpus`; local Grobid `http://localhost:8070`, version 0.8.1 |

Read-only status commands **on erenna**:

```bash
python -m json.tool /tmp/corpus-v15-full-candidate-20260922/runs/20260922T174832.798244Z/run.json
tail -n 40 /tmp/corpus-v15-full-candidate-20260922/runs/20260922T174832.798244Z/extract.log
ps -p 820851 -o pid,ppid,stat,etime,args
```

The driver runs normal CLI phases serially: taxonomy, extract, embed, post,
bundle. Every phase must pass its audit before the next starts. It checks a
50 GiB free-space floor before phases, holds an exclusive lock, records failures,
and preserves output. Do not start a competing extraction/embedding pipeline,
delete partial outputs, splice vectors, or bypass a failed audit.

If it stops, inspect the failure and verify no supervisor/children remain before
using the reviewed driver for normal implicit resume. The original launch was:

```bash
/home/claude/miniforge3/envs/corpus/bin/python -u \
  /tmp/corpus-v15-full-candidate-20260922/driver.py \
  --execute --authorization \
  /tmp/corpus-v15-full-candidate-20260922/launch-authorization.json
```

That is a recovery reference, **not a command to launch while this run is
active**. No renewed user permission is needed just to resume authorized work;
the driver authorization records reviewed input/code hashes and completed gold
acceptance. Changes to those inputs require a fresh, explicitly identified
candidate. The scratch README's "prepared, not launched" line is historical;
live `runs/` receipts supersede it.

Frozen inputs are independent, read-only copies of every retained reference
PDF, with full SHA-256 verification, plus BibTeX, lexicon, taxonomy and
instructions. CPU extraction uses four threads and OCR two jobs. The configured
local vision backend is deliberately deferred during extract-only processing
under #341. **This run supplies no Qwen acceptance.** Full launch identities and
hashes are committed in
[the launch receipt](examples/siphonophore_full_candidate_2026_09_22.json).

## Evidence already complete

- [CPU gold acceptance](examples/siphonophore_gold_refresh_2026_09_22/README.md):
  all 35 documents completed normal extraction, real BGE-M3 embedding, post,
  bundle, source scoring, live MCP and unchanged resume. All 3,336 vectors and
  chunk payloads were checked. Fourteen of 675 scored pages improve and the
  remaining page metrics are unchanged; source limitations remain documented.
  Figure/caption metrics are unchanged, not proof of correct Qwen crops.
- Unchanged gold extraction/embedding skip all documents. LanceDB adds a
  no-row deletion transaction during orphan pruning, so database file hashes
  differ; every row, vector, metadata field and generation is exactly equal.
  Acceptance is semantic equality, not byte-identical database files.
- [Current scientific notation review](../tests/fixtures/text_integrity/source_review/current_ab6353f/README.md)
  supports #303's closure. Keep its omitted/indeterminate source cases and
  selected-sample limitations; do not turn it into corpus-wide precision.
- [Named source pipeline pilots](examples/source_pages_2026_09_22/README.md)
  and [source acceptance inventory](V1_5_SOURCE_ACCEPTANCE.md) preserve the
  completed original issue evidence. Source fixes need not be rediscovered.
- Full local T0 previously passed: 2,767 tests, 26 optional skips, 90
  corpus/resume deselections. Required Ruff and subsequent focused checks pass.
  On `384911a`, all hosted PR checks passed: lint/unit, Linux, macOS, compose,
  and clean-room installation. Recheck CI for later commits; a handoff-only
  edit does not require repeating the full CPU gold build.

## Retrieval and bibliography: preserve these distinctions

Use [RETRIEVAL_EVALUATION.md](RETRIEVAL_EVALUATION.md), the exact
[frozen 22-query supplement](examples/siphonophore_retrieval_review_2026_09_22/README.md),
and the [guarded retained baseline](examples/siphonophore_retrieval_baseline_2026_09_22/README.md).
Twelve fixed queries have been reviewed. The ten independent queries' ranks,
distances, metrics and aggregate outcomes remain sealed until the candidate and
its controls are ready; do not inspect them opportunistically or change labels
after seeing rankings. The candidate policy was fixed before launch.

The retained full reference is **v1.2.1**, not the audited September 9
`1.4.0.dev0` bundle (pipeline SHA `734be4cfc301531bdcf90a1d23c2956bcbeed2fd`).
Its historic vector producer is unknown; matching the BGE-M3 model name and
dimension does not prove identical model weights. Comparison requires both
artifact and actual indexed-paper populations to match, with stable table
versions around queries. The retained zero-chunk paper `62e07c061591`
(Wangersky_Lane1960) is absent from its index. If the candidate recovers it,
report the membership difference; do not drop a competitor to pass the guard.

The [retained Pugh split receipt](examples/siphonophore_pugh_split_v1_2_1_2026_09_22.json)
reproduces the historical split, but its raw citation fields are empty. It
cannot justify bulk remapping or replace #296/#314's missing audited history.
Fresh named source cases are separate evidence with their own provenance.

## Assets outside Git and machine-specific constraints

All paths here are local to erenna. Before retiring that machine or clearing
`/tmp`, preserve needed evidence and outputs; a clone alone is insufficient.

| Asset | Local path |
| --- | --- |
| Source library and bibliography | `/home/claude/repos/siphonophores/library`; `/home/claude/repos/siphonophores/siphonophores.bib` |
| Library revision used | `5ad0164e6840ea16adb6e54a3a5712e20b234feb` |
| Completed gold, baseline, config and raw phase logs | `/tmp/corpus-v15-acceptance/gold` |
| Gold live-serving receipts | `/tmp/corpus-gold-serving-verification-20260922/run-2` |
| Retained v1.2.1 bundle | `/home/claude/corpuscles/siphonophore_20260831/_serve` |
| Guarded baseline raw captures, including sealed outcomes | `/tmp/corpus-issue320-retained-guarded-capture-20260922` |
| Full candidate, frozen inputs, scripts and live receipts | `/tmp/corpus-v15-full-candidate-20260922` |

The committed receipts identify local payload hashes; they do not contain every
large artifact or the scratch driver/audit scripts. Retrieve those from erenna
if continuing execution elsewhere, and preserve input identities and sealed
outcomes. The driver/audit contain absolute host paths: relocation requires
reviewing path changes and refreshing the corresponding authorization hashes.
Recreate the pinned code checkout from Git; copying a linked worktree alone
leaves its `.git` pointer referring to the old repository. Match the recorded
OCR/model dependencies and Grobid service before attempting a resumed build.
If relocating an unfinished build, stop it deliberately only when a transfer plan is ready;
do not treat a copy of a live mutable index as a completed acceptance bundle.

On this host, put the conda environment's `bin` directory first on `PATH`:
using its absolute Python executable alone is insufficient for OCR/test
subprocesses. LanceDB operations hang in the restricted sandbox and were run
with escalation. The GTX 1080 is unsupported by the installed Torch build;
there is no usable Bouchet SSH connection. Those are local limitations, not
product requirements. No fresh GPU pilot has been completed.

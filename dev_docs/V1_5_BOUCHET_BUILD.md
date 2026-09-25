# v1.5 candidate rebuild on Bouchet — September 25, 2026

This continues the [September 24 handoff](V1_5_HANDOFF.md). The user requested
a fresh build here, using SLURM for intensive work, and identified
`../siphonophores` as the source library. This is a new candidate; the erenna
run's state has not been established and its outputs have not been moved.

## September 25 recovery after the first status check

The `27479902` extraction chain was stopped after repeated runtime failures:
PyMuPDF `1.24.1` lacks `Pixmap.pil_image`, used during PDF preparation and
extraction. One batch also failed to start Grobid because its administrative
port was occupied. The original wrapper selected ports in Linux's ephemeral
range. The replacement uses the existing `slurm/bouchet_paths.sh` formula,
`8100 + 2 * (job_id % 400)`, and checks both ports on all interfaces.

All jobs in that chain were stopped before changing its wrapper. Its artifacts
remain in `output-attempt2/`, and its scripts and ledger are preserved with
`attempt2` names. A fresh `output/` contains the verified taxonomy copy. Frozen
PDFs and supporting inputs are unchanged.

Recovery job `27480182` installs **PyMuPDF 1.28.0**, matching erenna, into the
candidate's isolated `runtime/` directory and tests the required raster API.
It does not modify the shared conda environment. Two normal-pipeline source
pilots, array `27480190` (indices 65 and 261, batch size 1), replay the observed
extraction failure in AshaDevi_etal2010 and PDF-preparation failure in Hissmann2005.
Their completion requires clean summaries and an extraction completion record.
Launcher `27480191` submits the complete build chain only after both pass.
These recovery jobs were queued at this update; inspect their logs before
claiming recovery succeeded. `runtime-receipt.json` records the installed overlay;
the original `preparation.json` remains evidence of the initial environment.
The next `slurm-jobs.json` supersedes the cancelled job IDs below.

## Build identity and location

- Source library: `/nfs/roberts/project/pi_cwd7/cwd7/siphonophores`, revision
  `5ad0164e6840ea16adb6e54a3a5712e20b234feb` at preparation.
- Pinned code: `302bfac4d5e917424e607388fde3a04c0f3c1f00`, in
  `corpus/scratch/v15-bouchet-20260925/build-code`. The detached worktree keeps
  later development separate from running jobs.
- Candidate root:
  `/nfs/roberts/scratch/pi_cwd7/cwd7/siphonophore_v15_candidate_20260925`.
  Project storage had only about 175 GiB available; this scratch allocation
  had about 9.9 TiB available when selected.
- Preparation job: `27479572`. Automatic chain-submission job: `27479578`.
  An earlier pending preparation job, `27479569`, was cancelled before it ran
  to move the Grobid readiness check into each extraction allocation.

Preparation and chain submission completed successfully. The new candidate has
1,775 verified unique PDFs totaling 6,878,588,072 bytes, matching the frozen
September 22 candidate's population and source-byte count. These are dated
acceptance-run identities, not permanent library-size claims. Taxonomy ingestion
completed. The corrected extraction chain passed preflight and entered the
normal extraction module after the startup issue recorded below. The
[launch receipt](examples/siphonophore_bouchet_candidate_2026_09_25.json) pins
inputs, configuration, package/model versions, control-script hashes and jobs.

| Phase | SLURM job | Dependency |
| --- | --- | --- |
| Extraction | `27479902`, array `0–27%8` | preparation |
| Extraction audit | `27479904` | successful complete extraction array |
| GPU embedding and marker audit | `27479905` | extraction audit |
| Post, bundle, and SSE smoke | `27479906` | embedding |

The first chain (`27479864`–`27479867`) was stopped after task 0 passed
environment/Grobid preflight but the top-level CLI rejected `--strict-network`
before extraction. The corrected wrapper uses the already-exported, supported
`CORPUS_STRICT_NETWORK=1` setting. Its predecessor and job ledger remain as
`control/phase-attempt1.py` and `slurm-jobs-attempt1.json`. Frozen inputs and the
taxonomy database were retained; no source or production artifact was changed.

Docling, Torch, Transformers, sentence-transformers, LanceDB and PyArrow match
the erenna receipt's versions. Bouchet has PyMuPDF `1.24.1`, versus erenna's
`1.28.0`; this is recorded in `preparation.json` and must be considered in fresh
source validation. This launch does not establish equivalence of their output.

Preparation copies and verifies every retained reference PDF against its full
SHA-256, preserving the audited bundle's document population. It freezes separate
read-only copies of the bibliography, lexicon, taxonomy archive and instructions.
`preparation.json` records the exact PDF population, supporting-file hashes,
configuration, installed package versions, cached model revisions and language
packs. It must exist and preparation must finish successfully before extraction.

The configuration follows the completed September 22 gold configuration. CPU
extraction uses four model threads and two OCR workers. The accelerator setting
is `auto`: extraction allocations hide GPUs, while embedding requires an
allocated, supported GPU. Vision remains deferred, as in the earlier candidate.
The separate Qwen pilot below does not count as a full vision build.

The extraction array uses `ceil(unique PDFs / 64)` tasks, with at most eight
running simultaneously. Each task has its own Grobid service inside the same
allocation; readiness and strict network checks precede extraction. This avoids
a shared service expiring while late batches are still queued. Full extraction
completion, required stage records, document membership and recorded failures
are checked before embedding. Embedding completion markers gate post-processing
and bundling. Finalization checks bundle membership and runs the live SSE smoke
test. These operational checks do not replace scientific source acceptance or
the frozen retrieval evaluation.

## Monitor and continue

The preparation log is under
`corpus/scratch/v15-bouchet-20260925/build/prepare-27479572.{out,err}`.
After preparation, the candidate root contains copied control scripts, a frozen
configuration and `preparation.json`. The submission job records each downstream
job immediately in `slurm-jobs.json`; phase logs go into `logs/`, audit results
into `audits/`, and build artifacts into `output/`.

```bash
squeue -u cwd7
sacct -j 27479572,27479578 --format=JobID,State,ExitCode,Elapsed
cat /nfs/roberts/scratch/pi_cwd7/cwd7/siphonophore_v15_candidate_20260925/slurm-jobs.json
```

Inspect actual job receipts before resubmitting. Preparation refuses an existing
candidate directory; the chain submitter refuses an existing job ledger. Resume
failed extraction through normal implicit resume after checking that no writer
remains. Preserve completed output and the frozen inputs. Do not expose the ten
sealed independent retrieval outcomes until the candidate and its controls are
ready. The older retained v1.2.1 embedding-producer uncertainty still applies.

## Independent acceptance work recovered here

The local bundle at
`/nfs/roberts/project/pi_cwd7/cwd7/corpuscles/siphonophore_20260908_cea2b06/corpus_bundle`
matches the previously unavailable audited snapshot: `1.4.0.dev0`, created
`2026-09-09T09:23:59Z`, pipeline `734be4c`. Its original build artifacts and TEI
are also present outside the served bundle.

- **#337:** SLURM job `27478218` passed every SSE smoke-test layer against the
  separate `hydrozoa_20260920/corpus_bundle` production bundle, including actual
  query embedding. The corpus has no lexicon categories, so that optional check
  was explicitly skipped. This supplies the missing non-reference production
  transport check; it is not source-fidelity acceptance or a new fixture corpus.
- **#296/#314:** job `27479106` audits the historical bibliography and runs
  current authority/reconciliation against separate clean/incremental databases.
  These reuse historical extraction artifacts, not the fresh candidate.
  Historical Mańko authors are correct and ordered; the BibTeX stamp is absent,
  with no stored reconciliation decision for that paper. The reported Pugh
  populations reproduce 59 and 23 distinct citing papers with zero overlap.
  Unlike the v1.2.1 snapshot, this database retains raw observations. Inspect
  the job's final receipts before claiming repaired mappings or equivalence.
- **#305/#342:** H200 job `27479041` captured fresh model inputs, responses and
  crops for the three named cases and six additional panel controls. Five cases
  passed the tool's numeric checks; four did not, so the job exited 1. Scientific
  visual review remains incomplete. The reconstructed panel population contains
  1,868 figures rather than the original pilot's 1,865; this does not establish
  reproduction of that original six-case sample. Initial job `27478219` stopped
  in preparation after mistakenly including embedded figures in the sample;
  its partial directory is preserved separately.

Raw logs, scripts, historical inventory and pilot captures are under
`corpus/scratch/v15-bouchet-20260925`. No issue closure, release, production
replacement or retrieval improvement is established by these job submissions.

# Running the corpus pipeline on Bouchet (YCRC)

Operational notes for running the full ~2000-paper siphonophore corpus on Bouchet. See [OVERVIEW.md](OVERVIEW.md) for pipeline architecture and `dunnlab-hpc` skill for Dunn-lab-wide YCRC conventions.

## One-time setup

### 1. Repo + corpus on project storage

All paths in this guide and the `batch_*.sh` scripts derive from a single env var, `$BOUCHET_PROJECT`, defined in [slurm/bouchet_paths.sh](../slurm/bouchet_paths.sh). The default is `/nfs/roberts/project/pi_cwd7/cwd7` (PI Casey Dunn's lab project storage). To run from a different project root, edit that one line — every batch script sources it.

```bash
export BOUCHET_PROJECT=/nfs/roberts/project/pi_cwd7/cwd7
cd "$BOUCHET_PROJECT"
git clone git@github.com:caseywdunn/corpus.git
git clone git@github.com:dunnlab/siphonophores.git

# git-lfs is not on Bouchet's default PATH — load the module first.
# If `module load git-lfs` fails, search for the right name with
# `module spider lfs`, or install into the conda env as a fallback:
#   conda install -n corpus -c conda-forge git-lfs
module load git-lfs
cd siphonophores
git lfs install --local                   # wires LFS hooks into .git/config
git lfs pull                              # ~2000 PDFs, several GB
```

Resulting layout (all under `$BOUCHET_PROJECT`):

```
corpus/                       ← this repo (the `corpus` CLI)
siphonophores/                ← input PDF tree + siphonophores.bib + lexicon.yaml
siphonophore_corpuscle/       ← the corpuscle: config.yaml + build artifacts
  config.yaml                 ← authored in step 3 (the source of truth)
  documents/<HASH>/…          ← per-paper artifacts (created by extract)
  *.sqlite, vector_db/        ← cross-paper DBs + LanceDB (created by embed/post)
  _serve/                     ← distilled served bundle (created by bundle)
cache/huggingface/            ← model cache (see below)
```

As of #138 the SLURM phase scripts drive the **same `corpus run` CLI**
users run, one phase per job. All per-corpuscle inputs (PDF dir, BibTeX,
lexicon, taxonomy source, Grobid) live in `siphonophore_corpuscle/config.yaml`
— **not** as CLI flags or env vars. The scripts reference it through
`$CORPUS_CONFIG` (default `$BOUCHET_PROJECT/siphonophore_corpuscle/config.yaml`,
set in [slurm/bouchet_paths.sh](../slurm/bouchet_paths.sh)).

### 2. Conda environment

```bash
module load miniconda
conda env create -f "$BOUCHET_PROJECT/corpus/environment.yaml"
conda activate corpus
```

`environment.yaml` ships the Tesseract packs that back the default `ocr.ocr_languages_default` set in `config.yaml` (eng, deu, fra, rus, lat, spa, por, chi_sim, chi_tra, jpn, ell, kor, plus opt-in grc). The one exception is **19th-c. German Fraktur (`deu_latf`)**, which isn't on conda-forge — install it manually per [INSTALL.md](INSTALL.md#additional-ocr-language-packs).

**GPU torch (issue #21).** PyPI's default `torch` resolves to a CUDA 13 build, but Bouchet's H200 driver caps at CUDA 12.8. Re-install torch from PyTorch's cu128 index after the conda env exists:

```bash
conda activate corpus
pip install torch==2.9.0 torchvision==0.24.0 \
    --index-url https://download.pytorch.org/whl/cu128 \
    --force-reinstall
```

Verify with `python -c "import torch; print(torch.__version__, torch.cuda.is_available())"` on a `gpu` or `gpu_h200` node — should print `2.9.0+cu128 True`. The Stage 1 batch script's preflight (`slurm/batch_process_corpus.sh`) aborts loudly if torch / docling can't import, so a botched install fails fast instead of producing mis-structured output.

### 3. Author the corpuscle config.yaml

This is the source of truth for the build — `corpus run` reads every
per-corpuscle input from it. Scaffold it once, then edit:

```bash
conda activate corpus
mkdir -p "$BOUCHET_PROJECT/siphonophore_corpuscle"
cd "$BOUCHET_PROJECT/siphonophore_corpuscle"
corpus init                              # drops a commented config.yaml here
```

Edit `config.yaml` so it points at the siphonophores repo. Paths resolve
**relative to the config file's directory** (here, `siphonophore_corpuscle/`),
so the `../siphonophores/...` paths below work regardless of
`$BOUCHET_PROJECT` (YAML has no env-var expansion — keep them relative or
write absolute literals):

```yaml
# $BOUCHET_PROJECT/siphonophore_corpuscle/config.yaml
input_pdfs: ../siphonophores/library          # the PDF tree (was $INPUT_DIR)
output_dir: .                                 # artifacts land here in the
                                              #   corpuscle dir (= old $OUTPUT_DIR);
                                              #   keeps resume pointed at any
                                              #   existing documents/<HASH>/ tree
bib: ../siphonophores/siphonophores.bib       # was $BIB_FILE
lexicon: ../siphonophores/lexicon.yaml        # was $LEXICON

taxonomy:                                     # Siphonophorae, WoRMS AphiaID 1371
  source: worms
  root_id: 1371

figures:
  panel_detection: vision-local               # Pass 3b uses local Qwen2.5-VL

grobid:
  url: http://localhost:8070                  # overridden at submit time by
                                              #   $GROBID_URL (see step 5 / #138)
```

Validate it before submitting anything — `corpus check` confirms the host
can run the build and the config parses; `corpus run --dry-run` plans the
phases without writing artifacts:

```bash
corpus -c "$BOUCHET_PROJECT/siphonophore_corpuscle/config.yaml" check
corpus -c "$BOUCHET_PROJECT/siphonophore_corpuscle/config.yaml" run --dry-run
```

`$CORPUS_CONFIG` already defaults to this path in `bouchet_paths.sh`, so the
phase scripts find it with no extra flags. To build a *different* corpuscle,
export `CORPUS_CONFIG=/path/to/other/config.yaml` before submitting.

### 4. Pre-download HuggingFace models

Do this once on an interactive node (login nodes usually block outbound internet; request a compute node with `salloc -p interactive -t 0:30:00`):

```bash
# Pre-download into the SAME cache the batch jobs read —
# $BOUCHET_PROJECT/cache/huggingface, the path slurm/bouchet_paths.sh
# exports as HF_HOME. Don't use ~/project/... here: if it doesn't
# resolve to $BOUCHET_PROJECT, the jobs won't find these weights and
# re-download them (the "Stale HF_HOME" pitfall below).
export HF_HOME="$BOUCHET_PROJECT/cache/huggingface"
export TRANSFORMERS_CACHE="$HF_HOME/hub"
# Silence the implicit-token warning (#97). The batch jobs get this for
# free (it's set in pipeline/__init__.py on import), but these python -c
# one-liners import sentence_transformers/transformers directly and skip
# that, so set it here too.
export HF_HUB_DISABLE_IMPLICIT_TOKEN=1

# BGE-M3 (embeddings) — ~2 GB
python -c "from sentence_transformers import SentenceTransformer; \
    SentenceTransformer('BAAI/bge-m3')"

# Qwen2.5-VL-7B (Pass 3b local backend) — ~15 GB
python -c "from transformers import AutoProcessor, \
    Qwen2_5_VLForConditionalGeneration; \
    AutoProcessor.from_pretrained('Qwen/Qwen2.5-VL-7B-Instruct'); \
    Qwen2_5_VLForConditionalGeneration.from_pretrained('Qwen/Qwen2.5-VL-7B-Instruct')"
```

Batch jobs read the same `HF_HOME` (set by `slurm/bouchet_paths.sh` to
`$BOUCHET_PROJECT/cache/huggingface`), so they pick up these weights
without re-downloading.

**Is the pre-download idempotent?** *Partly.* Re-running these commands is
safe and a near-no-op when nothing changed upstream — `huggingface_hub`
HEAD-checks each file's etag against the cache and fetches only what
differs. But `from_pretrained(...)` / `SentenceTransformer(...)` resolve
the model's **`main` revision**, not a pinned commit, so the download is
**not version-locked**: if `BAAI/bge-m3` or the Qwen repo publishes a new
revision upstream, the next run pulls it into a new cache snapshot. That
can change behavior — a changed bge-m3 silently makes fresh query vectors
inconsistent with a `vector_db/lancedb` index built from the old weights
(same dimension, different values), and a changed Qwen alters Pass 3b
output.

Check what's cached and which commit each model currently resolves to:

```bash
huggingface-cli scan-cache                          # repos, REVISION (commit), refs, size, path
cat "$HF_HOME/hub/models--BAAI--bge-m3/refs/main"    # the commit `main` points at right now
```

For a reproducible build, **pin a commit** — pass `revision="<sha>"` to
`SentenceTransformer(..., revision=...)` / `from_pretrained(...,
revision=...)` and record it — and/or set `HF_HUB_OFFLINE=1` in the batch
environment so a job can never reach out and pull a newer revision: it
uses exactly the snapshot cached here or fails loudly.

### 5. Grobid as a Singularity service

Grobid is the extract phase's bibliographic + section-structure extractor. The `docker-compose.yml` in this repo targets macOS dev; on Bouchet run it via Singularity:

```bash
cd "$BOUCHET_PROJECT/cache"
singularity build grobid.sif docker://lfoppiano/grobid:0.8.1

# Start as a long-running job on an interactive node or via a separate
# SLURM job. Example interactive launch (adjust -t to cover stage 1):
salloc -p day -c 4 --mem=16G -t 24:00:00 \
    singularity run --bind "$BOUCHET_PROJECT" grobid.sif
# …returns a URL like http://<compute_node>:8070

# In a second terminal, submit the extract phase pointing at it. The
# extract script reads input/output/bib/lexicon from $CORPUS_CONFIG;
# $GROBID_URL overrides the config's grobid.url for this dynamically-
# allocated node (#138). Exported only when set, so config.yaml stays
# authoritative on standalone submits.
export GROBID_URL=http://<compute_node>:8070
sbatch slurm/batch_process_corpus.sh
```

In practice you don't submit Grobid by hand — `slurm/batch_pipeline.sh`
(see [Production run](#production-run)) starts it, discovers the node,
waits for `/api/isalive`, exports `$GROBID_URL` into the extract job, and
tears Grobid down afterward.

## Dry run (20–50 papers) before the full corpus

Before committing to the 2000-paper run, smoke-test end-to-end on a small
sample. Under the config-driven flow (#138) a sample is just a second
corpuscle — its own directory + `config.yaml` pointing at a slice of PDFs —
selected by exporting `CORPUS_CONFIG`:

```bash
# 1. A handful of PDFs spanning born-digital modern + historical scans +
#    German Fraktur to exercise the OCR / language / figure paths.
mkdir -p "$BOUCHET_PROJECT/siphonophores_sample/library"
# …copy ~30 PDFs into siphonophores_sample/library/

# 2. A sample corpuscle config. Copy the real one and repoint input_pdfs;
#    output_dir: . keeps the sample's artifacts out of the production tree.
mkdir -p "$BOUCHET_PROJECT/siphonophore_sample_corpuscle"
cd "$BOUCHET_PROJECT/siphonophore_sample_corpuscle"
cp "$BOUCHET_PROJECT/siphonophore_corpuscle/config.yaml" .
#   edit:  input_pdfs: ../siphonophores_sample/library
export CORPUS_CONFIG="$BOUCHET_PROJECT/siphonophore_sample_corpuscle/config.yaml"

# 3. Run the phases. The simplest path is the orchestrator, which inherits
#    $CORPUS_CONFIG via --export=ALL and runs the whole chain:
bash slurm/batch_pipeline.sh

#    …or submit phases by hand (each reads $CORPUS_CONFIG):
# sbatch slurm/batch_process_corpus.sh     # extract (CPU, day)
# sbatch slurm/batch_pass3b.sh             # vision  (GPU, gpu_h200)
# sbatch slurm/batch_embed.sh              # embed   (GPU)
# sbatch slurm/batch_finalize.sh           # post + bundle (CPU)
```

Inspect `$BOUCHET_PROJECT/siphonophore_sample_corpuscle/documents/<HASH>/summary.json` and `figures_report.html` for a handful of papers. Confirm:

- Grobid metadata populated (`metadata.json`)
- Figures extracted + captioned, panels listed in `figures.json`
- Pass 3b ROIs populated, compounds detected where expected
- Pass 3c renamed compound images to range notation (e.g., `fig_3-4.png`)
- `chunks.json` has reasonable content, embedded into `vector_db/lancedb`

Only then submit the three production jobs against the full input.

## Production run

The pipeline orchestrator (`slurm/batch_pipeline.sh`) handles Grobid startup, extract job-array submission, Grobid cleanup, and vision + embed + finalize chaining automatically — every job runs `corpus run --only <phase>` against `$CORPUS_CONFIG`.

The orchestrator parallelizes both Stage 1 (CPU) and Pass 3b (GPU)
independently. Stage 1 dominates wallclock for typical corpora; Pass 3b
is the long pole only when many figures need vision-LLM ROI detection.

| Knob | Stage | Default | What it controls |
|---|---|---|---|
| `NUM_BATCHES` | Stage 1 | `1` | Parallel CPU array tasks |
| `BATCH_SIZE` | Stage 1 | `256` | PDFs per Stage 1 task |
| `NUM_PASS3B_BATCHES` | Pass 3b | `1` | Parallel GPU array tasks |
| `PASS3B_BATCH_SIZE` | Pass 3b | `256` | Papers per Pass 3b task |

```bash
# Full pipeline with parallel Stage 1 batches of 256 PDFs each.
# Adjust NUM_BATCHES = ceil(total_unique_pdfs / 256).
NUM_BATCHES=8 bash slurm/batch_pipeline.sh

# Custom Stage 1 batch size:
NUM_BATCHES=4 BATCH_SIZE=512 bash slurm/batch_pipeline.sh

# Parallelize Pass 3b too — useful when figure-heavy corpora make the
# vision pass the long pole. NUM_PASS3B_BATCHES = ceil(total_papers /
# PASS3B_BATCH_SIZE). gpu_h200 partition has limited slots, so check
# `sinfo -p gpu_h200` before fanning out aggressively.
NUM_BATCHES=8 NUM_PASS3B_BATCHES=4 bash slurm/batch_pipeline.sh
```

The orchestrator:
1. Starts Grobid as a SLURM job, waits for it to be alive, exports `$GROBID_URL` for the extract job
2. Submits **extract** (`batch_process_corpus.sh`) as a job array (`--array=0-N`), each task processing a deterministic slice of the sorted hash list
3. Schedules Grobid cleanup after extract completes (runs regardless of success/failure)
4. Queues **vision** (`batch_pass3b.sh`, GPU array) and **embed** (`batch_embed.sh`, GPU) with `afterok` dependency on the full extract array
5. Queues **finalize** (`batch_finalize.sh` — `corpus run --only post` then `--only bundle`) with `afterok` dependency on embed

So the orchestrator now runs the whole build, including the cross-paper DBs and served-bundle distill, hands-off. Every job invokes `corpus -c "$CORPUS_CONFIG" run --only <phase>`.

For manual submission without the orchestrator (each phase reads `$CORPUS_CONFIG`; resume is implicit):

```bash
export GROBID_URL=http://<grobid_node>:8070      # extract needs Grobid
sbatch --array=0-7 slurm/batch_process_corpus.sh # extract
# After all array tasks complete:
sbatch --array=0-3 slurm/batch_pass3b.sh         # vision; or omit --array for single-task
sbatch slurm/batch_embed.sh                      # embed
sbatch slurm/batch_finalize.sh                   # post + bundle → _serve/
```

Resume is implicit in `corpus run`, so restarts are cheap — re-queuing a phase re-processes only papers whose inputs changed. Without `NUM_BATCHES` / `NUM_PASS3B_BATCHES` (or set to 1), the corresponding phase runs as a single job, which is usually right for small corpora.

## Post-pipeline cross-paper databases + served bundle

The extract/vision/embed phases produce per-paper artifacts under
`<output_dir>/documents/<HASH>/`. The **post** phase then runs the four
cross-paper builds in dependency order, and the **bundle** phase distills
the served bundle. Both are wrapped by `slurm/batch_finalize.sh`, which the
orchestrator chains automatically after embed — so normally you do nothing
here. The four builds (no longer separate top-level scripts; #138 folds
them into `corpus run --only post`) are:

1. **biblio authority** — deduplicated works + citation graph (`bib.authority`)
2. **taxon mentions** — cross-paper taxon-mention index (`pipeline.taxon_mentions`)
3. **in-text citations** — TEI body → `intext_citations.json` per paper (`pipeline.intext_citations`)
4. **reconcile** — merge ghost cited-references onto corpus papers (`bib.reconcile`)

To run the tail by hand (e.g. after fixing a build issue):

```bash
export CORPUS_CONFIG="$BOUCHET_PROJECT/siphonophore_corpuscle/config.yaml"

# Cross-paper DBs only. ENRICH_BHL=1 adds Biodiversity Heritage Library
# coverage for pre-DOI refs (slow, rate-limited ~1 req/s — many hours on
# the week partition); omit for a fast build without BHL.
sbatch slurm/batch_finalize.sh
# ENRICH_BHL=1 sbatch slurm/batch_finalize.sh

# Cross-paper DBs but skip the bundle distill:
# SKIP_BUNDLE=1 sbatch slurm/batch_finalize.sh
```

`batch_finalize.sh` runs `corpus run --only post` then `corpus run --only
bundle`. The post phase honors `--force-rebuild*` only when you ask; a plain
re-run is idempotent. Confirm the cross-paper SQLites landed:

```bash
ls -la "$BOUCHET_PROJECT/siphonophore_corpuscle"/{taxonomy,biblio_authority,taxon_mentions}.sqlite
```

The bundle phase distills into `<output_dir>/_serve/` — for the production
corpuscle that's `$BOUCHET_PROJECT/siphonophore_corpuscle/_serve/`. This
replaces the old standalone `package_for_serve.py` / `SERVE_BUNDLE_DIR`
step: `_serve/` **is** the deployable artifact (path-scrubbed + audited +
manifested). Re-distill any time with:

```bash
corpus -c "$CORPUS_CONFIG" run --only bundle
cat "$BOUCHET_PROJECT/siphonophore_corpuscle/_serve/bundle_manifest.json"
```

`_serve/` is what gets uploaded to S3 and consumed by the EC2 deploy — see [DEPLOY.md](DEPLOY.md) §6.

## Partition reference (from dev_docs/PLAN.md §7)

| Phase (`--only`) | Script | Partition | GPU? | Walltime |
|---|---|---|---|---|
| `extract` (OCR + docling + Grobid + Pass 2.5) | `slurm/batch_process_corpus.sh` | `day` | no | 24 h |
| `vision` (Pass 3b + 3c, Qwen2.5-VL-7B) | `slurm/batch_pass3b.sh` | `gpu_h200` | 1 | 24 h |
| `embed` (BGE-M3) | `slurm/batch_embed.sh` | `gpu` | 1 | 4 h |
| `post` + `bundle` (cross-paper DBs + served bundle) | `slurm/batch_finalize.sh` | `day` | no | 12 h |

Each script runs `corpus -c "$CORPUS_CONFIG" run --only <phase>`. Adjust walltimes if the corpus grows beyond ~2000 papers.

## Common pitfalls

- **Login nodes can't download models.** HF and other downloads must run on a compute node via `salloc`. Trying to pre-cache from a login node will hit outbound-firewall blocks.
- **Stale `HF_HOME`.** If a job re-downloads a model, `HF_HOME` isn't being honored — check that the export in the SLURM script points to a path you actually populated.
- **Grobid URL.** SLURM compute nodes can't talk to your laptop's `localhost:8070` — `$GROBID_URL` (which overrides the config's `grobid.url`) must resolve to a host visible from the job's node. If the Grobid node goes down mid-run, subsequent papers get placeholder metadata; a re-run's implicit resume won't retry them unless their inputs changed — force it with `corpus run --only extract --re-process-flagged <gate>` or by deleting the affected `metadata.json`.
- **Config not found / wrong corpuscle.** Every phase script reads `$CORPUS_CONFIG`. If a job dies with "no config.yaml" or builds the wrong tree, confirm `$CORPUS_CONFIG` points at the intended `config.yaml` (default is the production corpuscle; a leftover `export CORPUS_CONFIG=…sample…` from a smoke test will silently redirect a production submit).
- **LFS on extract.** If extract can't see the full PDFs (only LFS pointers), re-run `git lfs pull` in `$BOUCHET_PROJECT/siphonophores`.
- **GPU partition trap for the vision phase.** `rtx_5000_ada` nodes on the `gpu` partition carry an older NVIDIA driver (`nvidia-smi` reports CUDA 12.8 but torch 2.9.0+cu128 rejects it as "too old") and silently fall back to CPU, where Qwen2.5-VL-7B is unusable. Always submit the vision phase to `gpu_h200 --gpus=h200:1`. A pre-flight `python -c "import torch; assert torch.cuda.is_available()"` in `slurm/batch_pass3b.sh` aborts the job early if this is ever regressed.
- **lmod module cache flakiness.** `module load CUDA/12.6.0` occasionally errors with `CUDA/12.6.0.lua: Empty or non-existent file` even though `module spider CUDA` lists it. `slurm/batch_pass3b.sh` now skips the explicit CUDA module entirely — torch 2.9.0+cu128 ships bundled CUDA userspace libs in `site-packages/nvidia/` and works without it on gpu_h200. If another script hits this, retry with `module --ignore-cache load …`.

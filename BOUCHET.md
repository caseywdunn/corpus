# Running the corpus pipeline on Bouchet (YCRC)

Operational notes for running the full ~2000-paper siphonophore corpus on Bouchet. Pairs with PLAN.md §9 ("Bouchet prep") and the three SLURM batch scripts in this repo: `batch_process_corpus.sh`, `batch_pass3b.sh`, `batch_embed.sh`.

See also: `dunnlab-hpc` skill for Dunn-lab-wide YCRC conventions (partitions, storage layout, general SLURM).

## One-time setup

### 1. Repo + corpus on project storage

All paths in this guide and the `batch_*.sh` scripts derive from a single env var, `$BOUCHET_PROJECT`, defined in [bouchet_paths.sh](bouchet_paths.sh). The default is `/nfs/roberts/project/pi_cwd7/cwd7` (PI Casey Dunn's lab project storage). To run from a different project root, edit that one line — every batch script sources it.

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
corpus/             ← this repo
siphonophores/      ← input PDF tree
output/             ← created by stage 1
cache/huggingface/  ← model cache (see below)
```

### 2. Conda environment

```bash
module load miniconda
conda env create -f ~/project/corpus/environment.yaml
conda activate corpus
```

Additional Tesseract language packs (German Fraktur, French, Russian, Latin) still need to be installed — see `README.md` for specifics.

### 3. Pre-download HuggingFace models

Do this once on an interactive node (login nodes usually block outbound internet; request a compute node with `salloc -p interactive -t 0:30:00`):

```bash
export HF_HOME=~/project/cache/huggingface
export TRANSFORMERS_CACHE=$HF_HOME/hub

# BGE-M3 (embeddings) — ~2 GB
python -c "from sentence_transformers import SentenceTransformer; \
    SentenceTransformer('BAAI/bge-m3')"

# Qwen2.5-VL-7B (Pass 3b local backend) — ~15 GB
python -c "from transformers import AutoProcessor, \
    Qwen2_5_VLForConditionalGeneration; \
    AutoProcessor.from_pretrained('Qwen/Qwen2.5-VL-7B-Instruct'); \
    Qwen2_5_VLForConditionalGeneration.from_pretrained('Qwen/Qwen2.5-VL-7B-Instruct')"
```

Batch jobs set `HF_HOME` to the same path so they see the cached weights without re-downloading.

### 4. Grobid as a Singularity service

Grobid is stage 1's bibliographic + section-structure extractor. The `docker-compose.yml` in this repo targets macOS dev; on Bouchet run it via Singularity:

```bash
cd "$BOUCHET_PROJECT/cache"
singularity build grobid.sif docker://lfoppiano/grobid:0.8.1

# Start as a long-running job on an interactive node or via a separate
# SLURM job. Example interactive launch (adjust -t to cover stage 1):
salloc -p day -c 4 --mem=16G -t 24:00:00 \
    singularity run --bind "$BOUCHET_PROJECT" grobid.sif
# …returns a URL like http://<compute_node>:8070

# In a second terminal, submit stage 1 pointing at it:
export GROBID_URL=http://<compute_node>:8070
sbatch batch_process_corpus.sh
```

Alternative: run Grobid in the same job as stage 1 by adding `singularity exec ... &` before the `python process_corpus.py` call, with a loop that waits for `/api/isalive`. More brittle; the two-job approach above is cleaner.

## Dry run (20–50 papers) before the full corpus

Before committing to the 2000-paper run, smoke-test the pipeline end-to-end on a small sample:

```bash
# Pick a subset by copying a slice of the input tree
mkdir -p "$BOUCHET_PROJECT/siphonophores_sample"
# …copy ~30 PDFs spanning born-digital modern + historical scans +
#    German Fraktur to exercise the OCR / language / figure paths

# Override the defaults from bouchet_paths.sh via --export.
SAMPLE_IN="$BOUCHET_PROJECT/siphonophores_sample"
SAMPLE_OUT="$BOUCHET_PROJECT/output_sample"

# Stage 1 (CPU, day partition)
sbatch --job-name=corpus-sample-stage1 \
    --export=ALL,INPUT_DIR=$SAMPLE_IN,OUTPUT_DIR=$SAMPLE_OUT \
    batch_process_corpus.sh

# Pass 3b + 3c (GPU, ~30 min for 30 papers)
sbatch --job-name=corpus-sample-pass3b \
    --export=ALL,INPUT_DIR=$SAMPLE_IN,OUTPUT_DIR=$SAMPLE_OUT \
    batch_pass3b.sh

# Embeddings (GPU, <5 min for 30 papers)
sbatch --job-name=corpus-sample-embed \
    --export=ALL,OUTPUT_DIR=$SAMPLE_OUT \
    batch_embed.sh
```

Inspect `~/project/output_sample/documents/<HASH>/summary.json` and `figures_report.html` for a handful of papers. Confirm:

- Grobid metadata populated (`metadata.json`)
- Figures extracted + captioned, panels listed in `figures.json`
- Pass 3b ROIs populated, compounds detected where expected
- Pass 3c renamed compound images to range notation (e.g., `fig_3-4.png`)
- `chunks.json` has reasonable content, embedded into `vector_db/lancedb`

Only then submit the three production jobs against the full input.

## Production run

```bash
# Assume Grobid is up at $GROBID_URL; all three jobs can be chained with -d afterok:
J1=$(sbatch --parsable batch_process_corpus.sh)
J2=$(sbatch --parsable --dependency=afterok:$J1 batch_pass3b.sh)
J3=$(sbatch --parsable --dependency=afterok:$J2 batch_embed.sh)
echo "Chained: stage1=$J1, pass3b=$J2, embed=$J3"
```

`--resume` in every script means restarts are cheap — re-queuing picks up where the previous run left off.

## Partition reference (from PLAN.md §7)

| Stage | Script | Partition | GPU? | Walltime |
|---|---|---|---|---|
| Stage 1 (OCR + docling + Grobid + Pass 2.5) | `batch_process_corpus.sh` | `day` | no | 24 h |
| Pass 3b + 3c (Qwen2.5-VL-7B) | `batch_pass3b.sh` | `gpu_h200` | 1 | 24 h |
| Embeddings (BGE-M3) | `batch_embed.sh` | `gpu` | 1 | 4 h |

Adjust walltimes if the corpus grows beyond ~2000 papers.

## Common pitfalls

- **Login nodes can't download models.** HF and other downloads must run on a compute node via `salloc`. Trying to pre-cache from a login node will hit outbound-firewall blocks.
- **Stale `HF_HOME`.** If a job re-downloads a model, `HF_HOME` isn't being honored — check that the export in the SLURM script points to a path you actually populated.
- **Grobid URL.** SLURM compute nodes can't talk to your laptop's `localhost:8070` — `GROBID_URL` must resolve to a host visible from the job's node. If the Grobid node goes down mid-run, subsequent papers get placeholder metadata; `--resume` won't retry them without `--no-resume` or manually deleting `metadata.json`.
- **LFS on stage 1.** If stage 1 can't see the full PDFs (only LFS pointers), re-run `git lfs pull` in `~/project/siphonophores`.
- **GPU partition trap for Pass 3b.** `rtx_5000_ada` nodes on the `gpu` partition carry an older NVIDIA driver (`nvidia-smi` reports CUDA 12.8 but torch 2.9.0+cu128 rejects it as "too old") and silently fall back to CPU, where Qwen2.5-VL-7B is unusable. Always submit Pass 3b to `gpu_h200 --gpus=h200:1`. A pre-flight `python -c "import torch; assert torch.cuda.is_available()"` in `batch_pass3b.sh` aborts the job early if this is ever regressed.
- **lmod module cache flakiness.** `module load CUDA/12.6.0` occasionally errors with `CUDA/12.6.0.lua: Empty or non-existent file` even though `module spider CUDA` lists it. `batch_pass3b.sh` now skips the explicit CUDA module entirely — torch 2.9.0+cu128 ships bundled CUDA userspace libs in `site-packages/nvidia/` and works without it on gpu_h200. If another script hits this, retry with `module --ignore-cache load …`.

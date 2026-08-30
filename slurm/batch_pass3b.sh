#!/bin/bash
#SBATCH --job-name=corpus-pass3b
#SBATCH --partition=gpu_h200
#SBATCH --gpus=h200:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm-pass3b-%A_%a.out
#SBATCH --error=logs/slurm-pass3b-%A_%a.err
#
# Pass 3b + 3c: vision-model-driven panel detection + compound resolution.
# Uses the local Qwen2.5-VL-7B backend on GPU — no API cost, no network.
#
# Run AFTER Stage 1 (batch_process_corpus.sh) has completed.
#
# SLURM array support (#27): set --array=0-N when submitting and pass
# BATCH_SIZE in the env to slice the work across N+1 GPU jobs. Single-job
# mode (no --array) is unchanged. Example for 8 H200s in parallel:
#     BATCH_SIZE=256 sbatch --array=0-7 batch_pass3b.sh
#
# Partition: gpu_h200, for VRAM headroom and throughput on Qwen2.5-VL-7B
# — a preference, not a hard requirement. The whole GPU fleet runs
# driver 580.159.04 / CUDA 13.0, and the pinned torch 2.12.0 runs on
# every card type (verified 2026-08-01 on h200, b200,
# rtx_pro_6000_blackwell, rtx_5000_ada, a40, l40s). The preflight below
# is cheap and catches any future regression, so it stays.
#
# Measured runtime: 1 h 24 m as a single job over 1,769 papers
# (2026-08-04 build), against the 24 h wall.
#
# Only figures whose caption declares more than one panel reach the VLM
# (the `len(panels) <= 1: continue` filter in _pass3b_annotate_rois,
# pipeline/figure_passes.py) — on that build, 934 of 21,789 figure
# records, 4.3%. The GPU work is therefore small; the wallclock is
# dominated by walking every document directory and rewriting
# figures.json, so it scales with corpus size, not figure count.
#
# An earlier estimate here read "~20 figures/paper × 2000 papers ≈ 11–17
# hours"; it assumed every figure goes to the GPU and is what made Pass
# 3b look like the pipeline's long pole. It is not. The job still walks
# every document directory to rewrite figures.json, so wall time scales
# with corpus size, but the GPU work does not.
#
# Consequently --array is rarely worth it here: batch_pipeline.sh
# defaults NUM_PASS3B_BATCHES to 1. Fan out only for a genuinely
# figure-bound corpus, and mind the 16-GPU per-user cap on gpu_h200.
#
# Usage:
#     sbatch batch_pass3b.sh

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
# Pre-download the model to $HF_HOME (set by bouchet_paths.sh) before
# submitting:
#     corpus prefetch --include-vision
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"
echo "HuggingFace cache: $HF_HOME"

# ── Environment ──────────────────────────────────────────────────────
# torch ships bundled CUDA userspace libs (see site-packages/nvidia/),
# so no explicit CUDA module is needed on any GPU partition. We
# previously loaded CUDA/12.6.0 mirroring vial_scan, but
# that module file occasionally vanishes from the 2024a tree's lmod
# cache (observed job 8485894, 2026-04-16) — skipping it eliminates
# that failure mode.
# `module reset`, not `module purge` (#252). sbatch --export=ALL propagates
# LOADEDMODULES / _LMFILES_, so a batch job starts believing miniconda is
# loaded; purge unloads it and the modulefile's hook calls `conda`, which is
# a shell function that is not exported and does not exist in a
# non-interactive shell. Result: two alarming lines at the head of every
# job's stderr -- `conda: command not found` and a CondaError -- for jobs
# that then run fine. That noise cost real time once already, sitting above
# the actual cause of a failed Stage 1 launch (#251) and drawing the first
# hypothesis. Purge also does not do what it implies: StdEnv is sticky and
# survives it. `module reset` restores that same sticky default, matches
# YCRC's documented convention, and emits one informational line.
module reset
module load miniconda/24.7.1
eval "$(conda shell.bash hook)"
conda activate corpus

cd "$REPO_DIR"
mkdir -p logs

# ── Run ──────────────────────────────────────────────────────────────
echo "Starting Pass 3b (local VLM) at $(date)"
echo "GPU: $(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo "Config: $CORPUS_CONFIG"

# Sanity check: abort if torch can't see the GPU (driver mismatch etc.)
python -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 2)" || {
    echo "ERROR: torch.cuda.is_available() == False. Aborting before Pass 3b."
    exit 2
}

# Array-task batch slicing (#27). Activates only under sbatch --array.
# In single-job mode, no batch flags are passed and the run covers all
# papers.
BATCH_SIZE="${BATCH_SIZE:-256}"
BATCH_ARGS=()
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    BATCH_ARGS=(--batch-index "$SLURM_ARRAY_TASK_ID" --batch-size "$BATCH_SIZE")
    echo "Array task $SLURM_ARRAY_TASK_ID, batch size $BATCH_SIZE"
fi

# Vision phase only (Pass 3b + 3c). The `--only vision` orchestrator step
# already forces --refresh-vision / --no-grobid / --no-taxa internally;
# we only pin the backend to the local Qwen2.5-VL with --figure-panels
# vision-local. Resume is implicit. No Grobid needed (#138).
echo "corpus -c $CORPUS_CONFIG run --only vision --figure-panels vision-local ${BATCH_ARGS[*]}"

corpus -c "$CORPUS_CONFIG" run --only vision \
    --figure-panels vision-local \
    "${BATCH_ARGS[@]}"

echo "Pass 3b + 3c completed at $(date)"

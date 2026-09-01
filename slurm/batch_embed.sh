#!/bin/bash
#SBATCH --job-name=corpus-embed
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --constraint="gpu:a40|gpu:a5000"
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm-embed-%j.out
#SBATCH --error=logs/slurm-embed-%j.err
#
# Stage 2: embed chunks into LanceDB using the local BGE-M3 model.
# 48 GB on an A40 is far more than the 568M-parameter BGE-M3 needs at
# fp16; the type pin is about kernels, not capacity.
#
# CONSTRAIN THE GPU TYPE — a bare `--gpus=1` with no constraint is what
# broke the 2026-08-31 build. The `gpu` partition mixes four cards, and
# the pinned torch (2.12.0+cu130) ships kernels for
# sm_75/80/86/90/100/120:
#
#     a40             sm_86  ✓        l40s           sm_89  ✗
#     a5000           sm_86  ✓        rtx_5000_ada   sm_89  ✗
#
# Every arch is built as `code=sm_X` with no `code=compute_X` alongside
# it, so there is no embedded PTX for an unsupported card to JIT from.
# The gap is structural, not a slow path. Confirm with:
#   python -c "import torch; print(torch.__config__.show())" | tr ' ' '\n' | grep gencode
#
# Land on an sm_89 card and pipeline.accelerator correctly refuses it
# (#198) and falls back to CPU — so the job holds a GPU at 0% utilization
# while embedding at ~1/10th speed. On 2026-08-31 that ran 75 min, got
# through 181 of 1775 documents, and was cancelled by YCRC's job_defense
# GPU-utilization policy before it could finish.
#
# The `gpu:<type>` node features let one job accept either supported card,
# which schedules far sooner than pinning a single type — the a40 nodes
# are routinely half-drained while the 12 a5000 nodes sit open. Other
# supported types live in their own partitions: h200, b200,
# rtx_pro_6000_blackwell. See #270 for why the type constraint is a
# workaround: a torch bump that drops sm_86 reintroduces this silently.
#
# Run AFTER Stage 1 (batch_process_corpus.sh) has completed.
#
# Estimated runtime: full 2000-paper corpus embeds in well under an hour
# on a single GPU. The 4-hour wall is generous headroom.
#
# Usage:
#     sbatch batch_embed.sh

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
# Pre-download BGE-M3 to $HF_HOME (set by bouchet_paths.sh) before
# submitting:
#     corpus prefetch
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"
echo "HuggingFace cache: $HF_HOME"

# ── Environment ──────────────────────────────────────────────────────
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
module load miniconda CUDA
conda activate corpus

cd "$REPO_DIR"
mkdir -p logs

# ── Run ──────────────────────────────────────────────────────────────
echo "Starting embedding pass at $(date)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo "Config: $CORPUS_CONFIG"

# Embed phase only (BGE-M3 → LanceDB). Resume is implicit (#138).
corpus -c "$CORPUS_CONFIG" run --only embed || {
    EC=$?
    # Bus error (135) or segfault (139) during CUDA teardown after
    # successful embedding is a known issue on RTX 5000 Ada nodes with
    # older drivers. The embeddings are already written to disk — treat
    # this as success.
    if [ $EC -eq 135 ] || [ $EC -eq 139 ]; then
        echo "WARNING: Bus error during cleanup (exit $EC) — embeddings are complete"
    else
        echo "ERROR: corpus run --only embed failed with exit code $EC"
        exit $EC
    fi
}

echo "Embedding completed at $(date)"

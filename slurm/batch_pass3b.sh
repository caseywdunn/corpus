#!/bin/bash
#SBATCH --job-name=corpus-pass3b
#SBATCH --partition=gpu_h200
#SBATCH --gpus=h200:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm-pass3b-%j.out
#SBATCH --error=logs/slurm-pass3b-%j.err
#
# Pass 3b + 3c: vision-model-driven panel detection + compound resolution.
# Uses the local Qwen2.5-VL-7B backend on GPU — no API cost, no network.
#
# Run AFTER Stage 1 (batch_process_corpus.sh) has completed.
#
# Partition: we MUST use gpu_h200. The rtx_5000_ada nodes (gpu partition)
# carry an NVIDIA driver that torch 2.9.0+cu128 rejects as "too old"
# (reported driver 12080), causing a silent fallback to CPU that makes
# Qwen2.5-VL-7B unusable. H200 nodes (driver 570.x, supports CUDA 12.8)
# work out of the box. Confirmed 2026-04-16 via diag_gpu.sh on
# a1122u02n01 — all three module configs succeeded on H200.
#
# Estimated runtime: ~1–3 s/figure × ~20 figures/paper × 2000 papers ≈
# 11–17 hours. 24h wall gives headroom for the full corpus.
#
# Usage:
#     sbatch batch_pass3b.sh

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
# Pre-download the model to $HF_HOME (set by bouchet_paths.sh) on a
# compute node before submitting:
#     python -c "from transformers import AutoProcessor, Qwen2_5_VLForConditionalGeneration; \
#         AutoProcessor.from_pretrained('Qwen/Qwen2.5-VL-7B-Instruct'); \
#         Qwen2_5_VLForConditionalGeneration.from_pretrained('Qwen/Qwen2.5-VL-7B-Instruct')"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"
echo "HuggingFace cache: $HF_HOME"

# ── Environment ──────────────────────────────────────────────────────
# torch 2.9.0+cu128 ships bundled CUDA userspace libs (see
# site-packages/nvidia/), so no explicit CUDA module is needed on
# gpu_h200. We previously loaded CUDA/12.6.0 mirroring vial_scan, but
# that module file occasionally vanishes from the 2024a tree's lmod
# cache (observed job 8485894, 2026-04-16) — skipping it eliminates
# that failure mode.
module purge
module load miniconda/24.7.1
eval "$(conda shell.bash hook)"
conda activate corpus

cd "$REPO_DIR"
mkdir -p logs

# ── Run ──────────────────────────────────────────────────────────────
echo "Starting Pass 3b (local VLM) at $(date)"
echo "GPU: $(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo "Output: $OUTPUT_DIR"

# Sanity check: abort if torch can't see the GPU (driver mismatch etc.)
python -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 2)" || {
    echo "ERROR: torch.cuda.is_available() == False. Aborting before Pass 3b."
    exit 2
}

PY_ARGS=(
    "$INPUT_DIR" "$OUTPUT_DIR"
    --resume
    --refresh-vision
    --no-grobid
    --no-taxa
    --vision-backend local
)
echo "python args: process_corpus.py ${PY_ARGS[*]}"

python process_corpus.py "${PY_ARGS[@]}"

echo "Pass 3b + 3c completed at $(date)"

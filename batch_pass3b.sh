#!/bin/bash
#SBATCH --job-name=corpus-pass3b
#SBATCH --partition=gpu_h200
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm-pass3b-%j.out
#SBATCH --error=logs/slurm-pass3b-%j.err
#
# Pass 3b + 3c: vision-model-driven panel detection + compound resolution.
# Uses the local Qwen2.5-VL-7B backend on GPU — no API cost, no network.
#
# Run AFTER Stage 1 (batch_process_corpus.sh) has completed.
#
# The H200 partition is used because Qwen2.5-VL-7B at bfloat16 needs
# ~15 GB VRAM and benefits from the faster memory bandwidth. The RTX 5000
# (32 GB) on the gpu partition would also work — change --partition=gpu
# if H200 queues are long.
#
# Estimated runtime: ~1–3 s/figure × ~20 figures/paper × 2000 papers ≈
# 11–17 hours. The 12-hour wall is tight for the full corpus; increase
# to 24:00:00 if needed, or split the input into batches.
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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"
echo "HuggingFace cache: $HF_HOME"

# ── Environment ──────────────────────────────────────────────────────
module purge
module load miniconda CUDA
conda activate corpus

cd "$REPO_DIR"
mkdir -p logs

# ── Run ──────────────────────────────────────────────────────────────
echo "Starting Pass 3b (local VLM) at $(date)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo "Output: $OUTPUT_DIR"

python process_corpus.py \
    "$INPUT_DIR" "$OUTPUT_DIR" \
    --resume \
    --no-grobid \
    --no-taxa \
    --vision-backend local

echo "Pass 3b + 3c completed at $(date)"

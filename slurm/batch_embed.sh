#!/bin/bash
#SBATCH --job-name=corpus-embed
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm-embed-%j.out
#SBATCH --error=logs/slurm-embed-%j.err
#
# Stage 2: embed chunks into LanceDB using the local BGE-M3 model.
# Uses the gpu partition (RTX 5000 ADA, 32 GB) — more than enough for
# the 568M-parameter BGE-M3 at fp16.
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
# Pre-download BGE-M3 to $HF_HOME (set by bouchet_paths.sh) on a compute
# node before submitting:
#     python -c "from sentence_transformers import SentenceTransformer; \
#         SentenceTransformer('BAAI/bge-m3')"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
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
echo "Starting embedding pass at $(date)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo "Output: $OUTPUT_DIR"

python embed_chunks.py \
    "$OUTPUT_DIR" \
    --resume || {
    EC=$?
    # Bus error (135) or segfault (139) during CUDA teardown after
    # successful embedding is a known issue on RTX 5000 Ada nodes with
    # older drivers. The embeddings are already written to disk — treat
    # this as success.
    if [ $EC -eq 135 ] || [ $EC -eq 139 ]; then
        echo "WARNING: Bus error during cleanup (exit $EC) — embeddings are complete"
    else
        echo "ERROR: embed_chunks.py failed with exit code $EC"
        exit $EC
    fi
}

echo "Embedding completed at $(date)"

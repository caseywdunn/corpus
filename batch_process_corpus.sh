#!/bin/bash
#SBATCH --job-name=corpus-stage1
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm-stage1-%A_%a.out
#SBATCH --error=logs/slurm-stage1-%A_%a.err
#
# Stage 1: PDF processing pipeline (CPU-only).
# OCR, docling extraction, chunking, metadata (Grobid), Pass 2.5.
# No GPU needed — runs on the day partition.
#
# Usage:
#     # Ensure Grobid is running first (see batch_grobid.sh or Singularity notes)
#     sbatch batch_process_corpus.sh
#
# Paths are loaded from bouchet_paths.sh — edit that file to change
# $BOUCHET_PROJECT or any of the per-stage directories.
#
# Grobid: start as a Singularity service on an interactive node before
# submitting this job, or run the lightweight CRF image in a companion
# SLURM job:
#     singularity run --bind $BOUCHET_PROJECT docker://lfoppiano/grobid:0.8.1
# Then export GROBID_URL=http://<grobid-host>:8070

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
# sbatch copies this script to /var/spool/slurmd, so BASH_SOURCE can't
# locate bouchet_paths.sh. Use $SLURM_SUBMIT_DIR when running under
# sbatch; otherwise fall back to the script's own directory.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

# ── Environment ──────────────────────────────────────────────────────
module purge
module load miniconda
conda activate corpus

cd "$REPO_DIR"
mkdir -p logs "$OUTPUT_DIR"

# ── Grobid ───────────────────────────────────────────────────────────
GROBID_URL="${GROBID_URL:-http://localhost:8070}"
echo "Grobid URL: $GROBID_URL"

# Quick health check — warn but don't abort (pipeline has --no-grobid fallback)
if ! curl -fsS "$GROBID_URL/api/isalive" >/dev/null 2>&1; then
    echo "WARNING: Grobid not reachable at $GROBID_URL — metadata will use fallback"
fi

# ── WoRMS + anatomy lexicon (if available) ───────────────────────────
WORMS_FLAG=""
if [ -f "$OUTPUT_DIR/worms_siphonophora.sqlite" ]; then
    WORMS_FLAG="--worms-sqlite $OUTPUT_DIR/worms_siphonophora.sqlite"
fi
ANATOMY_FLAG=""
if [ -f "$REPO_DIR/resources/anatomy_lexicon.yaml" ]; then
    ANATOMY_FLAG="--anatomy-lexicon $REPO_DIR/resources/anatomy_lexicon.yaml"
fi

# ── Batch parameters (for SLURM job arrays) ─────────────────────────
BATCH_SIZE="${BATCH_SIZE:-256}"
BATCH_ARGS=""
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    BATCH_ARGS="--batch-index $SLURM_ARRAY_TASK_ID --batch-size $BATCH_SIZE"
    echo "Array task $SLURM_ARRAY_TASK_ID, batch size $BATCH_SIZE"
fi

# ── Run ──────────────────────────────────────────────────────────────
echo "Starting Stage 1 processing at $(date)"
echo "Input:  $INPUT_DIR"
echo "Output: $OUTPUT_DIR"

python process_corpus.py \
    "$INPUT_DIR" "$OUTPUT_DIR" \
    --resume \
    --grobid-url "$GROBID_URL" \
    $WORMS_FLAG \
    $ANATOMY_FLAG \
    $BATCH_ARGS

echo "Stage 1 completed at $(date)"

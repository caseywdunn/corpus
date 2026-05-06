#!/bin/bash
#SBATCH --job-name=corpus-stage1
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
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
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

# ── Environment ──────────────────────────────────────────────────────
module purge
module load miniconda
conda activate corpus

cd "$REPO_DIR"
mkdir -p logs "$OUTPUT_DIR"

# Fail-loud preflight (#20). docling import on Bouchet has historically
# fallen back to a broken tree-sitter / torch when LD_LIBRARY_PATH
# wasn't set right; the pipeline kept running and produced
# mis-structured output. Abort here instead.
echo "Preflight: import docling + torch ..."
python -c "
import sys
try:
    import docling.document_converter  # noqa: F401
    import torch  # noqa: F401
except Exception as e:
    sys.stderr.write(f'PREFLIGHT FAILED: {type(e).__name__}: {e}\n')
    sys.stderr.write(
        'Likely missing native libs. Confirm CONDA_PREFIX/lib is on '
        'LD_LIBRARY_PATH (handled by bouchet_paths.sh) and that '
        'torch matches the cluster CUDA version (issue #21).\n'
    )
    sys.exit(2)
print(f'docling + torch loaded OK ({torch.__version__})')
"

# ── Grobid ───────────────────────────────────────────────────────────
GROBID_URL="${GROBID_URL:-http://localhost:8070}"
echo "Grobid URL: $GROBID_URL"

# Quick health check — warn but don't abort (pipeline has --no-grobid fallback)
if ! curl -fsS "$GROBID_URL/api/isalive" >/dev/null 2>&1; then
    echo "WARNING: Grobid not reachable at $GROBID_URL — metadata will use fallback"
fi

# ── Taxonomy + lexicon (if available) ────────────────────────────────
# Taxonomy: built by ingest_taxonomy.py (any DwC source — WoRMS, GBIF,
# iNaturalist, or a curated CSV). Optional; pipeline degrades gracefully
# if absent. Lexicon: a single multi-category YAML; top-level keys are
# categories (anatomy, biogeography, …). See demo/lexicon.yaml.
TAXONOMY_FLAG=""
if [ -f "$OUTPUT_DIR/taxonomy.sqlite" ]; then
    TAXONOMY_FLAG="--taxonomy-db $OUTPUT_DIR/taxonomy.sqlite"
fi

LEXICON_FLAG=""
if [ -f "${LEXICON:-}" ]; then
    LEXICON_FLAG="--lexicon $LEXICON"
elif [ -f "$OUTPUT_DIR/lexicon.yaml" ]; then
    LEXICON_FLAG="--lexicon $OUTPUT_DIR/lexicon.yaml"
fi

BIB_FLAG=""
if [ -f "${BIB_FILE:-}" ]; then
    BIB_FLAG="--bib $BIB_FILE"
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
    $TAXONOMY_FLAG \
    $LEXICON_FLAG \
    $BIB_FLAG \
    $BATCH_ARGS

echo "Stage 1 completed at $(date)"

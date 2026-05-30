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
# $CORPUS_CONFIG is the source of truth for the Grobid address. When
# batch_pipeline.sh discovers the dynamically-allocated Grobid node it
# exports $GROBID_URL, which `corpus run` honors as an override (#138).
# Submitted standalone (no override), the grobid.url from config.yaml is
# used. Only export when explicitly set, so we never shadow the config.
if [ -n "${GROBID_URL:-}" ]; then
    export GROBID_URL
    echo "Grobid URL (override): $GROBID_URL"
    # Warn but don't abort — config may set grobid.disable, and the
    # pipeline degrades to header-fallback metadata if Grobid is down.
    if ! curl -fsS "$GROBID_URL/api/isalive" >/dev/null 2>&1; then
        echo "WARNING: Grobid not reachable at $GROBID_URL — metadata will use fallback"
    fi
else
    echo "Grobid URL: from $CORPUS_CONFIG (no \$GROBID_URL override)"
fi

# Inputs (taxonomy, lexicon, bib) are no longer passed as flags — they
# live in $CORPUS_CONFIG and `corpus run` reads them from there (#138).

# ── Batch parameters (for SLURM job arrays) ─────────────────────────
BATCH_SIZE="${BATCH_SIZE:-256}"
BATCH_ARGS=()
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    BATCH_ARGS=(--batch-index "$SLURM_ARRAY_TASK_ID" --batch-size "$BATCH_SIZE")
    echo "Array task $SLURM_ARRAY_TASK_ID, batch size $BATCH_SIZE"
fi

# ── Run ──────────────────────────────────────────────────────────────
# Extract phase only: OCR → docling → Grobid metadata → chunking → taxa +
# lexicon tagging. Resume is implicit; $GROBID_URL overrides the config's
# Grobid address (set above / by batch_pipeline.sh).
echo "Starting Stage 1 (extract) at $(date)"
echo "Config: $CORPUS_CONFIG"

corpus -c "$CORPUS_CONFIG" run --only extract "${BATCH_ARGS[@]}"

echo "Stage 1 completed at $(date)"

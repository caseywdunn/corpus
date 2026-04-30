#!/bin/bash
#SBATCH --job-name=biblio-authority
#SBATCH --partition=week
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/slurm-biblio-%j.out
#SBATCH --error=logs/slurm-biblio-%j.err
#
# Build the bibliographic authority database with BHL enrichment.
# Network + I/O bound (BHL API rate-limited at ~1 req/s with two-stage
# fallback), minimal CPU/memory. Uses the week partition for the long
# wall-clock time.
#
# Usage:
#     sbatch batch_biblio.sh
#
#     # Without BHL enrichment (fast, ~10 min):
#     BHL_ENRICH=0 sbatch batch_biblio.sh
#
#     # Incremental (skip --rebuild, only process new papers + retry BHL errors):
#     REBUILD=0 sbatch batch_biblio.sh

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

FINAL_DB="$OUTPUT_DIR/biblio_authority.sqlite"

# ── Environment ──────────────────────────────────────────────────────
module purge
module load miniconda
conda activate corpus

cd "$REPO_DIR"
mkdir -p logs

# ── Options ──────────────────────────────────────────────────────────
BHL_ENRICH="${BHL_ENRICH:-1}"
REBUILD="${REBUILD:-1}"

BHL_FLAG=""
if [ "$BHL_ENRICH" = "1" ]; then
    if [ -z "${BHL_API_KEY:-}" ]; then
        echo "ERROR: BHL_ENRICH=1 but BHL_API_KEY is not set."
        echo "  export BHL_API_KEY=<your-key>  # free at biodiversitylibrary.org/account"
        exit 1
    fi
    BHL_FLAG="--enrich-bhl --bhl-api-key $BHL_API_KEY"
fi

REBUILD_FLAG=""
if [ "$REBUILD" = "1" ]; then
    REBUILD_FLAG="--rebuild"
fi

# ── Run ──────────────────────────────────────────────────────────────
echo "Starting biblio authority build at $(date)"
echo "Output dir:  $OUTPUT_DIR"
echo "DB:          $FINAL_DB"
echo "BHL enrich:  $BHL_ENRICH"
echo "Rebuild:     $REBUILD"

python build_biblio_authority.py \
    "$OUTPUT_DIR" \
    -o "$FINAL_DB" \
    $REBUILD_FLAG \
    $BHL_FLAG \
    -v

echo "Biblio authority build completed at $(date)"

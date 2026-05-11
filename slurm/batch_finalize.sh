#!/bin/bash
#SBATCH --job-name=corpus-finalize
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm-finalize-%j.out
#SBATCH --error=logs/slurm-finalize-%j.err
#
# Finalize the cross-paper tail of a Bouchet build (#57).
#
# After Stage 1 + Stage 2 finish, four cross-paper builds still need to
# run in dependency order: build_biblio_authority → build_taxon_mentions
# → backfill_intext_citations → reconcile_corpus_to_biblio. v0.2 ran
# them on the login node manually after the array completed; v0.3
# wraps them in this single CPU job, which slurm/batch_pipeline.sh
# chains via --dependency=afterok: on the embed job.
#
# Single CPU job, no GPU. Mirrors batch_biblio.sh ergonomics + accepts
# the same opt-in env vars:
#
# Usage:
#     sbatch batch_finalize.sh
#
#     # With BHL enrichment (slow, network-bound):
#     ENRICH_BHL=1 sbatch batch_finalize.sh
#
#     # Distill into a served bundle when the build host ≠ serve host:
#     SERVE_BUNDLE_DIR=/path/to/bundle.v0.3.0 sbatch batch_finalize.sh

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

cd "$REPO_DIR"

ENRICH_BHL="${ENRICH_BHL:-0}"
SERVE_BUNDLE_DIR="${SERVE_BUNDLE_DIR:-}"

echo "=== Cross-paper finalize at $(date) ==="
echo "Output dir:  $OUTPUT_DIR"
echo "BHL enrich:  $ENRICH_BHL"
echo "Serve bundle: ${SERVE_BUNDLE_DIR:-<none>}"

BHL_FLAG=""
if [ "$ENRICH_BHL" = "1" ]; then
    BHL_FLAG="--enrich-bhl"
fi

# 1. Bibliographic authority DB (and apply any sidecar from #67)
echo "── build_biblio ──────────────────────────────────────────────"
python -m bib.authority "$OUTPUT_DIR" $BHL_FLAG -v

# 2. Taxon mention rollup
echo "── build_taxa ────────────────────────────────────────────────"
python -m pipeline.taxon_mentions "$OUTPUT_DIR" -v

# 3. In-text citation backfill from TEI bodies
echo "── backfill_intext ───────────────────────────────────────────"
python -m pipeline.intext_citations "$OUTPUT_DIR"

# 4. Reconcile ghost cited-references onto corpus papers
echo "── reconcile ─────────────────────────────────────────────────"
python -m bib.reconcile "$OUTPUT_DIR"

# 5. Optional: distill into a served bundle
if [ -n "$SERVE_BUNDLE_DIR" ]; then
    VERSION="$(python -c 'from pipeline.version import __version__; print(__version__)')"
    echo "── package_for_serve → $SERVE_BUNDLE_DIR (v$VERSION) ─────────"
    python -m mcpsrv.bundle "$OUTPUT_DIR" "$SERVE_BUNDLE_DIR" --version "$VERSION"
fi

echo "=== Cross-paper finalize complete at $(date) ==="

#!/bin/bash
#SBATCH --job-name=corpus-finalize
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
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
# Single CPU job, no GPU. Opt-in env vars:
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
module load miniconda
eval "$(conda shell.bash hook)"
conda activate corpus

ENRICH_BHL="${ENRICH_BHL:-0}"
SKIP_BUNDLE="${SKIP_BUNDLE:-0}"

echo "=== Cross-paper finalize at $(date) ==="
echo "Config:      $CORPUS_CONFIG"
echo "BHL enrich:  $ENRICH_BHL"
echo "Skip bundle: $SKIP_BUNDLE"

BHL_FLAG=()
if [ "$ENRICH_BHL" = "1" ]; then
    BHL_FLAG=(--enrich-bhl)
fi

# Post phase: the four cross-paper builds (biblio authority → taxon
# mentions → in-text citation backfill → reconcile) now run in dependency
# order inside `corpus run --only post` instead of four raw module
# invocations (#138).
echo "── post (cross-paper DBs) ────────────────────────────────────"
corpus -c "$CORPUS_CONFIG" run --only post "${BHL_FLAG[@]}"

# Bundle phase: distill the served bundle into <output_dir>/_serve/.
# This supersedes the old SERVE_BUNDLE_DIR / `mcpsrv.bundle` step — the
# served bundle now always lands at _serve/ (DEPLOY.md). Set SKIP_BUNDLE=1
# to stop after the cross-paper DBs.
if [ "$SKIP_BUNDLE" != "1" ]; then
    echo "── bundle (served-bundle distill → _serve/) ──────────────────"
    corpus -c "$CORPUS_CONFIG" run --only bundle
fi

echo "=== Cross-paper finalize complete at $(date) ==="

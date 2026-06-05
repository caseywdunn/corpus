#!/bin/bash
# launch_production.sh — verify sample corpuscle then kick off full production run.
# Submitted by batch_pipeline_watcher.sh with --dependency=afterok:<FINALIZE_JOB>
# so it only fires if the sample finalize completed successfully.
#
# Runs on the login node (via sbatch --partition=day, 5 min walltime).
#SBATCH --job-name=corpus-launch-prod
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/nfs/roberts/project/pi_cwd7/cwd7/corpus/logs/slurm-launch-prod-%j.out
#SBATCH --error=/nfs/roberts/project/pi_cwd7/cwd7/corpus/logs/slurm-launch-prod-%j.err

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

SAMPLE_CORPUSCLE="$BOUCHET_PROJECT/siphonophore_sample_corpuscle"
PROD_CORPUSCLE="$BOUCHET_PROJECT/siphonophore_corpuscle"

source /etc/profile.d/modules.sh
module load miniconda
eval "$(conda shell.bash hook)"
conda activate corpus

echo "=== Production launch check at $(date) ==="

# ── Verify sample output ──────────────────────────────────────────────
NDOCS=$(find "$SAMPLE_CORPUSCLE/documents" -name "summary.json" 2>/dev/null | wc -l)
echo "Sample documents with summary.json: $NDOCS"

if [ "$NDOCS" -lt 20 ]; then
    echo "ERROR: Only $NDOCS/30 sample papers have summary.json — sample run looks incomplete. Aborting production launch."
    exit 1
fi

HAS_VDB=0
[ -d "$SAMPLE_CORPUSCLE/vector_db/lancedb" ] && HAS_VDB=1
echo "Sample vector_db present: $HAS_VDB"

if [ "$HAS_VDB" -eq 0 ]; then
    echo "ERROR: No vector_db in sample corpuscle — embed phase may have failed. Aborting production launch."
    exit 1
fi

echo "Sample checks passed. Launching full production run."
echo ""

# ── Launch production pipeline ────────────────────────────────────────
export CORPUS_CONFIG="$PROD_CORPUSCLE/config.yaml"
cd "$REPO_DIR"

# 1763 PDFs / 256 per batch ≈ 7; use 8 for headroom.
# Check gpu_h200 availability; use 4 vision batches if slots exist.
H200_SLOTS=$(sinfo -p gpu_h200 -h -o "%A" 2>/dev/null | cut -d/ -f1 || echo 0)
echo "gpu_h200 available slots: $H200_SLOTS"

if [ "$H200_SLOTS" -ge 4 ]; then
    echo "Running: NUM_BATCHES=8 NUM_PASS3B_BATCHES=8 bash slurm/batch_pipeline.sh"
    NUM_BATCHES=8 NUM_PASS3B_BATCHES=8 bash slurm/batch_pipeline.sh
else
    echo "Running: NUM_BATCHES=8 NUM_PASS3B_BATCHES=1 bash slurm/batch_pipeline.sh (single vision batch)"
    NUM_BATCHES=8 NUM_PASS3B_BATCHES=1 bash slurm/batch_pipeline.sh
fi

echo "=== Production launch complete at $(date) ==="

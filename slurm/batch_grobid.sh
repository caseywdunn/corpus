#!/bin/bash
#SBATCH --job-name=grobid
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm-grobid-%j.out
#SBATCH --error=logs/slurm-grobid-%j.err
#
# Long-running Grobid service for stage 1 to talk to.
#
# Usage:
#     GROBID_JOB=$(sbatch --parsable batch_grobid.sh)
#     # wait for running, grab the node name, then export for stage 1:
#     NODE=$(squeue -j $GROBID_JOB -h -o %N)
#     export GROBID_URL=http://$NODE:8070
#     sbatch batch_process_corpus.sh
#
# Two binds matter:
#   1. $BOUCHET_PROJECT → so Grobid can read PDFs from the corpus tree.
#   2. A writable host dir → /opt/grobid/grobid-home/tmp. Without this
#      Grobid fails every request with HTTP 500 because the Singularity
#      rootfs is read-only and Grobid creates temp files per request.

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

cd "$REPO_DIR"
mkdir -p logs

# Per-job writable tmp so concurrent Grobid instances (if any) don't collide.
GROBID_TMP="$CACHE_DIR/grobid_tmp/$SLURM_JOB_ID"
mkdir -p "$GROBID_TMP"
trap 'rm -rf "$GROBID_TMP"' EXIT

echo "Grobid host: $(hostname)"
echo "Grobid tmp:  $GROBID_TMP"
echo "Starting at $(date)"

singularity run \
    --pwd /opt/grobid \
    --bind "$BOUCHET_PROJECT" \
    --bind "$GROBID_TMP:/opt/grobid/grobid-home/tmp" \
    "$CACHE_DIR/grobid.sif"

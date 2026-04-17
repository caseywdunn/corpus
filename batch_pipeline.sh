#!/bin/bash
# batch_pipeline.sh — orchestrate the full corpus processing pipeline.
#
# Chains: Grobid server → Stage 1 → cancel Grobid → Pass 3b + Embed
#
# This script runs on the login node (lightweight — just sbatch + curl +
# sleep). It auto-discovers the Grobid node, waits for the HTTP service,
# and wires up SLURM dependencies so the whole pipeline runs hands-off.
#
# Usage:
#     bash batch_pipeline.sh
#
#     # Override input/output for the sample dataset:
#     INPUT_DIR=$BOUCHET_PROJECT/siphonophores_sample \
#     OUTPUT_DIR=$BOUCHET_PROJECT/output_sample \
#         bash batch_pipeline.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

mkdir -p "$REPO_DIR/logs"

echo "=== Corpus Pipeline Launcher ==="
echo "Input:  $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo ""

# ── Step 1: Start Grobid server ─────────────────────────────────────
echo "Submitting Grobid server..."
GROBID_JOB=$(sbatch --parsable "$SCRIPT_DIR/batch_grobid.sh")
echo "  Grobid job: $GROBID_JOB"

# ── Step 2: Wait for Grobid to reach RUNNING state ──────────────────
echo "Waiting for Grobid to start..."
MAX_WAIT=600  # 10 minutes
ELAPSED=0
STATE="UNKNOWN"
while [ "$ELAPSED" -lt "$MAX_WAIT" ]; do
    STATE=$(squeue -j "$GROBID_JOB" -h -o "%T" 2>/dev/null || echo "UNKNOWN")
    if [ "$STATE" = "RUNNING" ]; then
        break
    fi
    sleep 10
    ELAPSED=$((ELAPSED + 10))
    printf "\r  Waited %ds (state: %s)..." "$ELAPSED" "$STATE"
done
echo ""

if [ "$STATE" != "RUNNING" ]; then
    echo "ERROR: Grobid job $GROBID_JOB did not start within ${MAX_WAIT}s (state: $STATE)"
    echo "       Cancel it with: scancel $GROBID_JOB"
    exit 1
fi

GROBID_NODE=$(squeue -j "$GROBID_JOB" -h -o "%N")
GROBID_URL="http://${GROBID_NODE}:8070"
echo "  Grobid running on $GROBID_NODE"

# ── Step 3: Wait for Grobid HTTP service to be ready ────────────────
echo "Waiting for Grobid HTTP service at $GROBID_URL ..."
GROBID_READY=0
for i in $(seq 1 60); do
    if curl -fsS "${GROBID_URL}/api/isalive" >/dev/null 2>&1; then
        GROBID_READY=1
        echo "  Grobid is alive (after ~$((i * 5))s)"
        break
    fi
    sleep 5
done

if [ "$GROBID_READY" -eq 0 ]; then
    echo "WARNING: Grobid not responding at $GROBID_URL after 5 minutes."
    echo "         Stage 1 will proceed but may fall back to no-grobid mode."
fi

# ── Step 4: Submit Stage 1 with Grobid URL ──────────────────────────
echo ""
BATCH_SIZE="${BATCH_SIZE:-256}"
NUM_BATCHES="${NUM_BATCHES:-1}"
echo "Submitting Stage 1 ($NUM_BATCHES batch(es) of $BATCH_SIZE)..."
if [ "$NUM_BATCHES" -gt 1 ]; then
    STAGE1_JOB=$(GROBID_URL="$GROBID_URL" sbatch --parsable \
        --array="0-$((NUM_BATCHES - 1))" \
        --export="ALL,BATCH_SIZE=$BATCH_SIZE" \
        "$SCRIPT_DIR/batch_process_corpus.sh")
else
    STAGE1_JOB=$(GROBID_URL="$GROBID_URL" sbatch --parsable \
        "$SCRIPT_DIR/batch_process_corpus.sh")
fi
echo "  Stage 1 job: $STAGE1_JOB"

# ── Step 5: Submit Grobid cleanup job (runs after Stage 1) ──────────
# afterany: runs whether stage1 succeeds or fails, so Grobid is always
# cleaned up. Uses a minimal allocation.
CANCEL_JOB=$(sbatch --parsable \
    --dependency=afterany:"$STAGE1_JOB" \
    --partition=day --time=00:05:00 --cpus-per-task=1 --mem=1G \
    --job-name=grobid-cancel \
    --output="$REPO_DIR/logs/slurm-grobid-cancel-%j.out" \
    --wrap="scancel $GROBID_JOB 2>/dev/null && echo 'Cancelled Grobid job $GROBID_JOB' || echo 'Grobid job $GROBID_JOB already finished'")
echo "  Grobid cancel job: $CANCEL_JOB (runs after Stage 1)"

# ── Step 6: Submit Pass 3b and Embed (depend on Stage 1) ────────────
echo ""
echo "Submitting Pass 3b (vision, GPU)..."
PASS3B_JOB=$(sbatch --parsable --dependency=afterok:"$STAGE1_JOB" "$SCRIPT_DIR/batch_pass3b.sh")
echo "  Pass 3b job: $PASS3B_JOB"

echo "Submitting Embed (GPU)..."
EMBED_JOB=$(sbatch --parsable --dependency=afterok:"$STAGE1_JOB" "$SCRIPT_DIR/batch_embed.sh")
echo "  Embed job: $EMBED_JOB"

# ── Summary ─────────────────────────────────────────────────────────
echo ""
echo "=== Pipeline Submitted ==="
echo ""
printf "  %-20s %s\n" "Grobid server:" "$GROBID_JOB (running on $GROBID_NODE)"
printf "  %-20s %s\n" "Stage 1:" "$STAGE1_JOB (running now)"
printf "  %-20s %s\n" "Grobid cancel:" "$CANCEL_JOB (after Stage 1)"
printf "  %-20s %s\n" "Pass 3b:" "$PASS3B_JOB (after Stage 1)"
printf "  %-20s %s\n" "Embed:" "$EMBED_JOB (after Stage 1)"
echo ""
echo "Monitor: squeue --me"
echo "Cancel all: scancel $GROBID_JOB $STAGE1_JOB $CANCEL_JOB $PASS3B_JOB $EMBED_JOB"

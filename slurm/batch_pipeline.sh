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
#     # Build a corpuscle other than corpuscles/current (e.g. the smoke test):
#     CORPUS_CONFIG=$BOUCHET_PROJECT/corpuscles/siphonophore_gold_YYYYMMDD/config.yaml \
#         bash batch_pipeline.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

mkdir -p "$REPO_DIR/logs"

echo "=== Corpus Pipeline Launcher ==="
echo "Config: $CORPUS_CONFIG"
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
BATCH_SIZE="${BATCH_SIZE:-64}"

# Auto-size the array from the corpuscle's own input_pdfs. The old
# NUM_BATCHES default of 1 meant a bare `bash batch_pipeline.sh` ran the
# whole corpus as a single task, and an operator-supplied value that was
# too small silently left the tail of the library unprocessed — the
# failure mode BOUCHET.md could only warn about in prose.
#
# Counting *files* over-provisions slightly relative to unique hashes
# (duplicate PDFs collapse), which is harmless: an empty batch exits 0 in
# a few minutes. Under-provisioning is the dangerous direction.
if [ -z "${NUM_BATCHES:-}" ]; then
    _cfg_dir="$(dirname "$CORPUS_CONFIG")"
    # Strip the key, any trailing ` # comment`, surrounding quotes and
    # trailing whitespace. Deliberately a grep and not a YAML parse: this
    # runs on the login node before any conda env is activated.
    _in_rel="$(sed -n 's/^input_pdfs:[[:space:]]*//p' "$CORPUS_CONFIG" | head -1 \
        | sed -e 's/[[:space:]]*#.*$//' -e 's/^["'"'"']//' -e 's/["'"'"']$//' \
              -e 's/[[:space:]]*$//')"
    if [ -z "$_in_rel" ]; then
        echo "ERROR: could not read input_pdfs from $CORPUS_CONFIG." >&2
        echo "       Set NUM_BATCHES explicitly to override." >&2
        exit 1
    fi
    _in_abs="$(cd "$_cfg_dir" && cd "$_in_rel" 2>/dev/null && pwd)" || {
        echo "ERROR: input_pdfs '$_in_rel' (from $CORPUS_CONFIG) is not a directory." >&2
        exit 1
    }
    _n_pdfs=$(find "$_in_abs" -iname '*.pdf' | wc -l)
    if [ "$_n_pdfs" -eq 0 ]; then
        echo "ERROR: no PDFs found under $_in_abs." >&2
        exit 1
    fi
    NUM_BATCHES=$(( (_n_pdfs + BATCH_SIZE - 1) / BATCH_SIZE ))
    echo "Auto-sized Stage 1: $_n_pdfs PDFs / $BATCH_SIZE per task = $NUM_BATCHES tasks"
fi

# Pre-build the taxonomy before anything is submitted (#251).
#
# The array tasks run `corpus run --only extract`, and that path does NOT
# build taxonomy.sqlite the way a full `corpus run` does — the orchestrator
# hard-errors before any work starts when it is missing. Two consecutive
# siphonophore builds lost their whole chain to this: `corpus check` had
# called the missing sqlite self-resolving, 26 of 28 array tasks died ~60s
# in, and `afterok` took Pass 3b, Embed and Finalize down with them.
#
# Doing it here, on the login node, rather than trusting the operator to
# have read a warning.
#
# Safe to repeat, which was asserted before it was true (#262): `names` had
# no uniqueness constraint and plain INSERTs, so this call doubled that table
# on every launch -- 801 rows to 1,602 to 2,403. Since 1.2.2 the ingest holds
# a unique index and INSERT OR IGNORE, and it deduplicates on open, so an
# unconditional call is genuinely a no-op *and* repairs a corpuscle damaged
# by 1.2.1 the next time this runs. About a second either way.
if [ -n "$CORPUS_CONFIG" ]; then
    echo "Ensuring taxonomy.sqlite is built..."
    if ! corpus taxonomy ingest; then
        echo "ERROR: taxonomy ingest failed. Stage 1 would fail the same way" >&2
        echo "       ~60s into every array task, taking the chain with it." >&2
        echo "       Fix the taxonomy config before resubmitting." >&2
        exit 1
    fi
fi

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
# Deliberately NOT defaulted to $NUM_BATCHES. Pass 3b only sends a figure
# to the VLM when its caption declares multiple panels
# (pipeline/figure_passes.py), which on the siphonophore corpus is 934 of
# 21,789 figure records — under 5%. A single GPU task covered the whole
# 1,769-paper library in 1 h 24 m. Fanning out to match Stage 1's array
# (now tens of tasks) would queue against the 16-GPU per-user cap on
# gpu_h200 for no gain. Raise it only for corpora that really are
# figure-bound.
NUM_PASS3B_BATCHES="${NUM_PASS3B_BATCHES:-1}"
PASS3B_BATCH_SIZE="${PASS3B_BATCH_SIZE:-256}"
echo "Submitting Pass 3b (vision, GPU; $NUM_PASS3B_BATCHES batch(es) of $PASS3B_BATCH_SIZE)..."
if [ "$NUM_PASS3B_BATCHES" -gt 1 ]; then
    PASS3B_JOB=$(sbatch --parsable \
        --dependency=afterok:"$STAGE1_JOB" \
        --array="0-$((NUM_PASS3B_BATCHES - 1))" \
        --export="ALL,BATCH_SIZE=$PASS3B_BATCH_SIZE" \
        "$SCRIPT_DIR/batch_pass3b.sh")
else
    PASS3B_JOB=$(sbatch --parsable \
        --dependency=afterok:"$STAGE1_JOB" \
        "$SCRIPT_DIR/batch_pass3b.sh")
fi
echo "  Pass 3b job: $PASS3B_JOB"

echo "Submitting Embed (GPU)..."
EMBED_JOB=$(sbatch --parsable --dependency=afterok:"$STAGE1_JOB" "$SCRIPT_DIR/batch_embed.sh")
echo "  Embed job: $EMBED_JOB"

# Cross-paper finalize (#57): the four post-pipeline builds plus the
# served-bundle distill.
#
# Depends on Pass 3b as well as Embed. Those two are siblings, so gating
# only on Embed let `bundle` start while Pass 3b was still rewriting
# figures.json and Pass 3c was still renaming split-panel PNGs — both of
# which mcpsrv/bundle.py copies into _serve/. The bundle could therefore
# capture pre-vision ROIs and stale figure filenames.
echo "Submitting Finalize (cross-paper tail)..."
FINALIZE_JOB=$(sbatch --parsable \
    --dependency=afterok:"$EMBED_JOB":"$PASS3B_JOB" \
    "$SCRIPT_DIR/batch_finalize.sh")
echo "  Finalize job: $FINALIZE_JOB"

# ── Step 7: Chain watchdog ──────────────────────────────────────────
# The 2026-08-02 build failed silently: two Stage 1 array tasks hit the
# wall, `afterok` was therefore never satisfied, and SLURM cancelled
# Pass 3b, Embed and Finalize with no log, no mail and no queue entry.
# afterany guarantees this job runs no matter how Stage 1 ends, so there
# is always one place that says what happened.
WATCHDOG_JOB=$(sbatch --parsable \
    --dependency=afterany:"$STAGE1_JOB" \
    --partition=day --time=00:10:00 --cpus-per-task=1 --mem=2G \
    --job-name=chain-watchdog \
    --output="$REPO_DIR/logs/slurm-watchdog-%j.out" \
    --wrap="set -u
echo '=== Stage 1 array outcome ==='
sacct -j $STAGE1_JOB --format=JobID%20,JobName%18,State%12,Elapsed,End%20
echo
n_done=\$(ls -1 '$(dirname "$CORPUS_CONFIG")'/documents/*/summary.json 2>/dev/null | wc -l)
echo \"Documents complete: \$n_done\"
if sacct -j $STAGE1_JOB -n -o State%20 | grep -qvE 'COMPLETED|^[[:space:]]*\$'; then
    echo
    echo '*** WARNING: at least one Stage 1 task did not COMPLETE. ***'
    echo '*** afterok was not satisfied, so Pass 3b ($PASS3B_JOB), Embed ($EMBED_JOB)'
    echo '*** and Finalize ($FINALIZE_JOB) have been CANCELLED by SLURM.'
    echo '*** Resume is implicit — resubmit with: bash slurm/batch_pipeline.sh'
else
    echo 'Stage 1 fully COMPLETED; downstream chain is live.'
fi")
echo "  Chain watchdog: $WATCHDOG_JOB (after Stage 1, always runs)"

# ── Summary ─────────────────────────────────────────────────────────
echo ""
echo "=== Pipeline Submitted ==="
echo ""
printf "  %-20s %s\n" "Grobid server:" "$GROBID_JOB (running on $GROBID_NODE)"
printf "  %-20s %s\n" "Stage 1:" "$STAGE1_JOB (running now)"
printf "  %-20s %s\n" "Grobid cancel:" "$CANCEL_JOB (after Stage 1)"
printf "  %-20s %s\n" "Pass 3b:" "$PASS3B_JOB (after Stage 1)"
printf "  %-20s %s\n" "Embed:" "$EMBED_JOB (after Stage 1)"
printf "  %-20s %s\n" "Finalize:" "$FINALIZE_JOB (after Embed + Pass 3b; #57)"
printf "  %-20s %s\n" "Chain watchdog:" "$WATCHDOG_JOB (after Stage 1, always)"
echo ""

# Show the dependency state so the chain is visibly live at submit time.
# Every downstream job should appear with a Dependency reason; a job
# missing from this list was rejected outright.
# `|| true` throughout: this block is diagnostic, and `set -e` must not
# abort a successfully-submitted pipeline over a transient squeue error.
echo "=== Dependency check ==="
squeue --me -o "%.12i %.18j %.10T %.30E" || true
echo ""
for _j in "$PASS3B_JOB" "$EMBED_JOB" "$FINALIZE_JOB"; do
    if ! squeue -j "$_j" -h -o "%i" >/dev/null 2>&1; then
        echo "WARNING: job $_j is not in the queue — it was never accepted." >&2
    fi
done

echo "Monitor: squeue --me"
echo "Watchdog log: $REPO_DIR/logs/slurm-watchdog-$WATCHDOG_JOB.out"
echo "Cancel all: scancel $GROBID_JOB $STAGE1_JOB $CANCEL_JOB $PASS3B_JOB $EMBED_JOB $FINALIZE_JOB $WATCHDOG_JOB"

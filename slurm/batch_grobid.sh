#!/bin/bash
#SBATCH --job-name=grobid
#SBATCH --partition=week
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/slurm-grobid-%j.out
#SBATCH --error=logs/slurm-grobid-%j.err
#
# Long-running Grobid service for stage 1 to talk to.
#
# `week`/48 h rather than `day`/24 h: this job is submitted *before*
# stage 1 and `day` caps at 24 h, so an equal wall guaranteed Grobid
# died first. Papers stage 1 processes after that get placeholder
# metadata, and implicit resume will NOT retry them — their inputs are
# unchanged — so the loss is silent. Outliving stage 1 costs nothing:
# the `afterany` grobid-cancel job in batch_pipeline.sh tears this down
# as soon as stage 1 ends, however it ends.
#
# Usage:
#     GROBID_JOB=$(sbatch --parsable batch_grobid.sh)
#     # wait for running, grab the node name, then export for stage 1:
#     until [ "$(squeue -j "$GROBID_JOB" -h -o %T)" = RUNNING ]; do sleep 5; done
#     source bouchet_paths.sh   # for corpus_grobid_port
#     PORT=$(corpus_grobid_port "$GROBID_JOB")
#     export GROBID_URL="http://$(squeue -j "$GROBID_JOB" -h -o %N):$PORT"
#     # RUNNING only means SLURM started the container. Grobid needs another
#     # ~30-60 s to load its models and bind the port, so poll before using it —
#     # a "Connection refused" in that window is startup, not failure.
#     # The port is derived from the job ID, not fixed at 8070 (#279); the job's
#     # own stdout echoes it as "Grobid URL:" if you would rather read it off.
#     until curl -fsS "$GROBID_URL/api/isalive" >/dev/null 2>&1; do sleep 5; done
#     sbatch batch_process_corpus.sh
#
# batch_pipeline.sh does all of the above for you; the manual path is for
# debugging and for the `corpus check` preflight (see dev_docs/BOUCHET.md §6).
#
# Three binds matter:
#   1. $BOUCHET_PROJECT → so Grobid can read PDFs from the corpus tree.
#   2. A writable host dir → /opt/grobid/grobid-home/tmp. Without this
#      Grobid fails every request with HTTP 500 because the Singularity
#      rootfs is read-only and Grobid creates temp files per request.
#   3. A writable host dir → /opt/grobid/logs. Without this Grobid crashes
#      on startup trying to open logs/grobid-service.log (FileNotFoundException).

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$SCRIPT_DIR/bouchet_paths.sh" ] || SCRIPT_DIR="$SCRIPT_DIR/slurm"
# shellcheck source=bouchet_paths.sh
source "$SCRIPT_DIR/bouchet_paths.sh"

cd "$REPO_DIR"
mkdir -p logs

# Per-job writable dirs so concurrent Grobid instances don't collide.
# Two binds required (comment in script header) plus a third:
#   3. A writable host dir → /opt/grobid/logs. Grobid tries to open
#      logs/grobid-service.log relative to --pwd and crashes when the
#      Singularity rootfs is read-only and the directory doesn't exist.
GROBID_TMP="$CACHE_DIR/grobid_tmp/$SLURM_JOB_ID"
GROBID_LOGS="$CACHE_DIR/grobid_logs/$SLURM_JOB_ID"

# Reap dirs leaked by earlier jobs. The trap below empties its dirs but often
# cannot remove them (see there), so sweep on the way in as well as on the way
# out. -mindepth 1 keeps the parents, which may not exist yet on a fresh
# project root.
#
# Age is the only safe signal for grobid_tmp: Grobid writes there per request
# and leaves it empty when idle, so a live instance's tmp dir looks exactly
# like a leaked one. --time=2-00:00:00 caps a live instance's age at two days,
# so -mtime +3 cannot touch a running job. Keep this threshold strictly above
# the walltime above if you ever change it.
find "$CACHE_DIR/grobid_tmp" "$CACHE_DIR/grobid_logs" -mindepth 1 -maxdepth 1 \
    -type d -mtime +3 -exec rm -rf {} + 2>/dev/null || true

# grobid_logs additionally admits a prompt rule, which is what actually keeps
# the leak in check: a running instance opens grobid-service.log within seconds
# of mkdir, so an empty log dir 10 min old is always dead.
find "$CACHE_DIR/grobid_logs" -mindepth 1 -maxdepth 1 \
    -type d -empty -mmin +10 -exec rmdir {} + 2>/dev/null || true

mkdir -p "$GROBID_TMP" "$GROBID_LOGS"

# Grobid still holds grobid-service.log open when SIGTERM lands, so NFS
# silly-renames it to .nfsXXXX and the rmdir hits "Directory not empty" until
# the client releases the file — measured at over 20 s, longer than SLURM's
# KillWait=30 s leaves us. So expect to lose that race and leave an empty dir
# behind; the sweep above reaps it on the next submit. What matters here is
# that the contents go and that cleanup never decides the job's exit status —
# the old `trap 'rm -rf ...'` failed loudly and returned 1 into a
# `set -e` EXIT trap.
cleanup() {
    rm -rf "$GROBID_TMP" "$GROBID_LOGS" 2>/dev/null || true
}
trap cleanup EXIT

echo "Grobid host: $(hostname)"
echo "Grobid tmp:  $GROBID_TMP"
echo "Starting at $(date)"

# A port of this job's own, so SLURM co-scheduling two of these on one node
# is no longer fatal (#279). Derived from the job ID by the shared function in
# bouchet_paths.sh; batch_pipeline.sh derives the same pair from the same job
# ID, so the client cannot drift from the server.
GROBID_PORT=$(corpus_grobid_port "${SLURM_JOB_ID:-0}")
GROBID_ADMIN_PORT=$(corpus_grobid_admin_port "${SLURM_JOB_ID:-0}")

# Both connectors have to move, and this is the part that is easy to get
# wrong: Dropwizard binds an *admin* connector besides the application one
# (default 8071), so overriding only the application port still dies with the
# same BindException. Verified against lfoppiano/grobid:0.8.1 — with just
# `applicationConnectors[0].port` changed, a second instance on a shared
# network exits 1 on `java.net.BindException: Address already in use`; with
# both changed, three instances ran side by side and returned byte-identical
# TEI for the same PDF.
#
# `dw.`-prefixed system properties are Dropwizard's own config-override
# mechanism, so this needs no change to Grobid or to the image. GROBID_SERVICE_OPTS
# in the image carries the java.library.path and --add-opens flags the service
# needs; JAVA_OPTS is separate and unset by default, so setting it clobbers
# nothing. SINGULARITYENV_ is the prefix both Singularity 3.x and Apptainer
# honour, unlike `--env`, whose support depends on the runtime version.
export SINGULARITYENV_JAVA_OPTS="-Ddw.server.applicationConnectors[0].port=$GROBID_PORT -Ddw.server.adminConnectors[0].port=$GROBID_ADMIN_PORT"

# Backstop for the derivation, not the main defence: two concurrent jobs land
# on the same pair only if their IDs are congruent mod 400. Say so plainly if
# it ever happens, because the downstream failure mode is nasty — this job has
# already reached RUNNING, so a pipeline waiting on job state alone would point
# Stage 1 at this node and be served by the *other* chain's Grobid.
if command -v ss >/dev/null 2>&1 && \
   ss -ltn 2>/dev/null | grep -qE ":($GROBID_PORT|$GROBID_ADMIN_PORT)[[:space:]]"; then
    echo "ERROR: port $GROBID_PORT or $GROBID_ADMIN_PORT is already bound on $(hostname)." >&2
    echo "       Another Grobid job derived the same pair from its job ID (#279)." >&2
    echo "       Resubmit — the next job ID will map elsewhere — or share a single" >&2
    echo "       server by passing GROBID_URL to batch_process_corpus.sh." >&2
    exit 1
fi

echo "Grobid port: $GROBID_PORT (admin $GROBID_ADMIN_PORT)"
echo "Grobid URL:  http://$(hostname):$GROBID_PORT"

singularity run \
    --pwd /opt/grobid \
    --bind "$BOUCHET_PROJECT" \
    --bind "$GROBID_TMP:/opt/grobid/grobid-home/tmp" \
    --bind "$GROBID_LOGS:/opt/grobid/logs" \
    "$CACHE_DIR/grobid.sif"

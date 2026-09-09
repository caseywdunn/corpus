# Single source of truth for Bouchet paths.
# Sourced by every batch_*.sh script in this repo, and referenced by
# BOUCHET.md for one-time setup commands.
#
# Nothing here names a specific corpuscle: switching builds must never
# require editing a tracked file. Point at a different build either by
# moving the `corpuscles/current` symlink or by exporting $CORPUS_CONFIG
# for one-off runs, e.g.:
#     CORPUS_CONFIG=$BOUCHET_PROJECT/corpuscles/siphonophore_gold_20260825/config.yaml \
#         bash batch_pipeline.sh

# Lab project storage (PI Casey Dunn). All persistent state lives under
# this path on Bouchet — never use $HOME for corpus data, since $HOME has
# a 125 GB quota that the model cache + LanceDB index will exceed.
BOUCHET_PROJECT="${BOUCHET_PROJECT:-/nfs/roberts/project/pi_cwd7/cwd7}"

REPO_DIR="${REPO_DIR:-$BOUCHET_PROJECT/corpus}"

# A console entry point uses its installed package location, not the shell's
# cwd. The orchestrator also changes cwd before starting Python subprocesses.
# Pin imports for both to the selected checkout, including alternate worktrees.
REPO_DIR="$(cd "$REPO_DIR" && pwd -P)" || return 1
export REPO_DIR
export PYTHONPATH="$REPO_DIR${PYTHONPATH:+:$PYTHONPATH}"

corpus_check_checkout() {
    # Run outside the repository root so cwd cannot mask a stale installation.
    # Call after activating conda, before each build phase or any submission.
    (
        cd "$REPO_DIR/pipeline" || exit 1
        python - "$REPO_DIR" "${CORPUS_BUILD_GIT_SHA:-}" <<'PY'
import importlib.util
from pathlib import Path
import subprocess
import sys

root = Path(sys.argv[1]).resolve()
for name in ("pipeline", "bib", "mcpsrv"):
    spec = importlib.util.find_spec(name)
    actual = Path(spec.origin).resolve() if spec and spec.origin else None
    expected = root / name / "__init__.py"
    if actual != expected:
        sys.exit(f"Checkout mismatch: {name} resolves to {actual}; expected {expected}")
    print(f"Build import: {name} = {actual}", flush=True)
sha = subprocess.check_output(
    ["git", "-C", str(root), "rev-parse", "HEAD"], text=True,
).strip()
if sys.argv[2] and sha != sys.argv[2]:
    sys.exit(f"Checkout changed since submission: {sha}; expected {sys.argv[2]}")
print(f"Build checkout: {root} @ {sha}", flush=True)
PY
    )
}

# Corpuscle config.yaml — the source of truth for the `corpus run` flow
# (#138). input_pdfs / output_dir / bib / lexicon / taxonomy / grobid all
# live INSIDE that file, not as CLI flags or shell vars here. The phase
# scripts (batch_process_corpus, batch_pass3b, batch_embed,
# batch_finalize) drive `corpus -c "$CORPUS_CONFIG" run --only <phase>`.
# See BOUCHET.md for how to author it. The dynamic Grobid node URL is
# still injected at submit time via $GROBID_URL, which `corpus run`
# honors as an override.
#
# Resolution order — deliberately no dated default, so starting a new
# build never means editing this tracked file:
#
#   1. $CORPUS_CONFIG, if exported (one-off runs, sample builds).
#      `corpus` itself honors this var too (cli.py `_resolve_config_path`,
#      #61), and sbatch propagates it by default (--export=ALL).
#   2. $BOUCHET_PROJECT/corpuscles/current — an untracked symlink naming
#      the active build. Retarget it to switch builds:
#          ln -sfn siphonophore_YYYYMMDD "$BOUCHET_PROJECT/corpuscles/current"
#      `ls -l "$BOUCHET_PROJECT/corpuscles/current"` then shows, on the
#      cluster rather than in git, which corpuscle jobs will write to.
#   3. Otherwise: hard error. Silently falling back to a stale corpuscle
#      would resume someone else's finished build.
if [[ -z "${CORPUS_CONFIG:-}" ]]; then
    _corpus_current="$BOUCHET_PROJECT/corpuscles/current/config.yaml"
    if [[ -e "$_corpus_current" ]]; then
        CORPUS_CONFIG="$_corpus_current"
    else
        echo "ERROR: no corpuscle selected." >&2
        echo "  \$CORPUS_CONFIG is unset and $BOUCHET_PROJECT/corpuscles/current" >&2
        echo "  does not resolve to a config.yaml. Pick a build with:" >&2
        echo "      ln -sfn siphonophore_YYYYMMDD \"$BOUCHET_PROJECT/corpuscles/current\"" >&2
        echo "  or export CORPUS_CONFIG=/path/to/corpuscle/config.yaml for a one-off." >&2
        unset _corpus_current
        # `return` when sourced (the normal case); `exit` if run directly.
        return 1 2>/dev/null || exit 1
    fi
    unset _corpus_current
fi
# Exported so batch_pipeline.sh pins one corpuscle for the whole chain:
# sbatch propagates it (--export=ALL), so every dependent job builds the
# corpuscle that was selected at submit time even if `current` is
# retargeted while Stage 1 is still running.
export CORPUS_CONFIG

CACHE_DIR="${CACHE_DIR:-$BOUCHET_PROJECT/cache}"

# HuggingFace cache (docling models + BGE-M3 + Qwen2.5-VL weights).
# Warm it once with `corpus prefetch` before the first submission — from
# the login node or anywhere else convenient. See BOUCHET.md.
#
# HF_HOME is the only cache var that matters. We deliberately do NOT set
# TRANSFORMERS_CACHE: transformers 5.x ignores it, so setting it only
# creates the illusion of a second knob.
export HF_HOME="${HF_HOME:-$CACHE_DIR/huggingface}"

# Native libs for the conda env (#20). Bouchet compute nodes don't
# always have $CONDA_PREFIX/lib on the runtime loader path; without
# it docling/torch silently fall back to whatever the system has and
# produce mis-structured output without crashing. Prepend the env's
# libs so they win — matches the resolution order pip / conda used
# at install time. Caller must `conda activate corpus` before sourcing.
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
fi

# ── Per-job Grobid ports (#279) ──────────────────────────────────────
# Grobid's Dropwizard service binds a fixed 8070, and SLURM is free to
# co-schedule several of these jobs on one node — at which point every
# instance after the first dies ~10 s in with a Jetty BindException.
# Submitting six pipelines put five Grobid jobs on two nodes and three
# failed. Worse than a plain failure: the job has already reached RUNNING,
# so a chain that waits on job state alone points Stage 1 at that node and
# is served by *another chain's* server. That happened, and two documents
# recorded `grobid.outcome = extracted` from the wrong service.
#
# So derive a port pair from the job ID. Two rules, both load-bearing:
#
#   * **Both connectors must move.** Dropwizard binds an *admin* connector
#     as well, default 8071, and overriding only the application port
#     still dies on it — verified against lfoppiano/grobid:0.8.1, which is
#     the trap the issue's own suggested fix would have hit. Stride of 2
#     with admin = app + 1 means no two jobs' pairs can overlap.
#   * **Stay clear of 8070/8071**, so a hand-started default-port Grobid —
#     the documented debugging path — never collides with a batch job.
#
# 400 slots. Two concurrent jobs collide only if their IDs are congruent
# mod 400, and concurrently submitted IDs differ by small amounts; the
# `ss` preflight in batch_grobid.sh is the backstop either way.
#
# Callers must derive rather than hardcode, and batch_pipeline.sh derives
# from the same job ID as the server, so the two cannot drift.
corpus_grobid_port() {
    local job_id="${1:-0}"
    echo $(( 8100 + 2 * (job_id % 400) ))
}

corpus_grobid_admin_port() {
    echo $(( $(corpus_grobid_port "${1:-0}") + 1 ))
}

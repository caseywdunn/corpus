#!/usr/bin/env bash
# dev_docs/ec2_smoke.sh
#
# One-shot platform-portability smoke test on a clean Ubuntu EC2 host.
# Validates the linux-x86_64 install path end-to-end: apt deps,
# miniforge, conda env from environment.yaml, pip install -e .,
# Grobid via Docker, the demo `corpus run`, bundle distillation, and
# the SSE/MCP round-trip. Exits 0 iff every success criterion in
# dev_docs/PLATFORM_SMOKE.md passes.
#
# Pre-flight (do these yourself; the script picks up from a corpus
# clone in cwd):
#
#   1. Launch an EC2 instance.
#        AMI:           Ubuntu 24.04 LTS (HVM, EBS-backed) — Canonical's
#                       official "ubuntu-noble-24.04-amd64-server-…" image
#        Instance type: c6i.2xlarge (8 vCPU, 16 GB) — fastest extract for
#                       ~$0.34/hr. m6i.xlarge works at ~$0.19/hr if you
#                       don't mind the slower extract phase.
#        EBS:           80 GB gp3 (Grobid full image is ~32 GB; we use
#                       the CRF-only image at ~7 GB here, so headroom
#                       is comfortable — the size is for safety, not
#                       necessity)
#        Security group: SSH (22) inbound from your IP only
#
#   2. SSH in.
#
#   3. Clone the repo and cd into it (uncomment to run, or do it by
#      hand — the script doesn't run these for you so it stays
#      idempotent against re-runs from inside an existing clone):
#
#        # sudo apt-get update -qq && sudo apt-get install -y -qq git
#        # git clone https://github.com/caseywdunn/corpus.git
#        # cd corpus
#
#   4. Run this script:
#
#        bash dev_docs/ec2_smoke.sh
#
# The script writes:
#   ./demo/output/                — pipeline artifacts (gitignored)
#   /tmp/corpus_ec2_smoke.log     — full `corpus run` output
#   /tmp/corpus_ec2_smoke.summary — machine-readable pass/fail tallies
#
# Teardown: terminate the EC2 instance when done. Single-use, nothing
# to preserve.
#
# Total wall time: ~20–25 min (Grobid pull + corpus run dominate). The
# WoRMS taxonomy ingest in `corpus run` was narrowed from Siphonophorae
# (order, ~10 min) to Cystonectae (suborder, ~2 min) in commit 7d0bd98
# specifically to keep this script under ~25 min.

set -eo pipefail
# NOTE: `set -u` would break conda's activate script (conda references
# unset variables internally). Keep nounset off.

# ── Sanity: are we in a clone of the corpus repo? ─────────────────
if [ ! -f environment.yaml ] || [ ! -f pyproject.toml ] \
   || ! grep -q 'name = "corpus"' pyproject.toml; then
    echo "ERROR: run this script from the root of a corpus repo clone." >&2
    echo "Expected: ./environment.yaml, ./pyproject.toml with name=\"corpus\"." >&2
    exit 2
fi

REPO_ROOT="$(pwd)"
LOG=/tmp/corpus_ec2_smoke.log
SUMMARY=/tmp/corpus_ec2_smoke.summary

# ── Reporting helpers ─────────────────────────────────────────────
if [ -t 1 ]; then
    GREEN=$'\033[1;32m'; RED=$'\033[1;31m'; YELLOW=$'\033[1;33m'
    CYAN=$'\033[1;36m'; RESET=$'\033[0m'
else
    GREEN=''; RED=''; YELLOW=''; CYAN=''; RESET=''
fi
PASS_COUNT=0
FAIL_COUNT=0
FAIL_DETAILS=()

section() { printf '\n%s═══ %s ═══%s\n' "$CYAN" "$*" "$RESET"; }
note_pass() {
    printf '  %s[PASS]%s %s\n' "$GREEN" "$RESET" "$1"
    PASS_COUNT=$((PASS_COUNT + 1))
}
note_fail() {
    printf '  %s[FAIL]%s %s\n' "$RED" "$RESET" "$1"
    FAIL_COUNT=$((FAIL_COUNT + 1))
    FAIL_DETAILS+=("$1")
}
note_warn() { printf '  %s[warn]%s %s\n' "$YELLOW" "$RESET" "$1"; }

phase_started_at=$(date +%s)
elapsed() {
    local now diff
    now=$(date +%s); diff=$((now - phase_started_at))
    printf '%s(phase took %dm %02ds)%s\n' "$YELLOW" \
        "$((diff / 60))" "$((diff % 60))" "$RESET"
    phase_started_at=$now
}

# ── Phase 1: system tools + Docker ────────────────────────────────
section "Phase 1 — apt deps + Docker"
sudo apt-get update -qq

# Required: things the pipeline genuinely can't run without.
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    git curl ca-certificates build-essential \
    ghostscript tesseract-ocr pandoc jq

# Optional ocrmypdf helpers — runtime degrades gracefully when these
# are missing (pipeline/scan.py drops --optimize 2→1 when pngquant
# isn't on PATH; INSTALL.md documents both as system-install extras).
# Some Ubuntu AMIs don't have `universe` enabled by default, so try
# best-effort and warn rather than aborting the whole smoke for a
# helper we already designed around.
for opt in pngquant jbig2enc; do
    if sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq "$opt" 2>/dev/null; then
        :
    else
        note_warn "$opt not available via apt — output PDFs will be larger, but smoke still runs"
    fi
done

# Docker via upstream apt repo (more current than Ubuntu's docker.io).
# Skip if already installed (re-run safety).
if ! command -v docker >/dev/null; then
    sudo install -m 0755 -d /etc/apt/keyrings
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg \
        | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
    sudo chmod a+r /etc/apt/keyrings/docker.gpg
    . /etc/os-release
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
         https://download.docker.com/linux/ubuntu ${VERSION_CODENAME} stable" \
        | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
    sudo apt-get update -qq
    sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        docker-ce docker-ce-cli containerd.io
fi
note_pass "system tools + Docker installed"
elapsed

# ── Phase 2: miniforge ────────────────────────────────────────────
section "Phase 2 — miniforge (conda)"
if [ ! -d "$HOME/miniforge3" ]; then
    curl -fsSL -o /tmp/Miniforge3.sh \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash /tmp/Miniforge3.sh -b -p "$HOME/miniforge3"
    rm /tmp/Miniforge3.sh
fi
# Source conda for this script; `conda init` only writes to .bashrc,
# which non-interactive shells don't pick up.
# shellcheck source=/dev/null
source "$HOME/miniforge3/etc/profile.d/conda.sh"
note_pass "miniforge ready ($(conda --version))"
elapsed

# ── Phase 3: corpus conda env + editable install ──────────────────
section "Phase 3 — conda env + pip install -e ."
if ! conda env list | awk '{print $1}' | grep -qx corpus; then
    conda env create -f environment.yaml
else
    note_warn "conda env 'corpus' already exists — reusing (script is idempotent)"
fi
conda activate corpus
pip install --quiet -e .

# Architecture sanity (would catch x86_64-under-Rosetta on macOS;
# on Ubuntu it's just a label).
host_arch=$(python -c "import platform; print(platform.machine())")
note_pass "corpus env active; host arch = $host_arch"
elapsed

# ── Phase 4: tessdata + pyflakes lint precheck ────────────────────
section "Phase 4 — tessdata + pyflakes precheck"
bash tools/install_tessdata.sh > /tmp/tessdata.out 2>&1 \
    && note_pass "tessdata language packs installed" \
    || { note_fail "tessdata install failed"; tail -20 /tmp/tessdata.out; }

# Fast-fail signal that the source tree compiles cleanly before we
# spend 15+ min on the demo run.
if python -m pytest tests/test_no_undefined_names.py -q \
       > /tmp/pyflakes.out 2>&1; then
    note_pass "pyflakes gate (tests/test_no_undefined_names.py)"
else
    note_fail "pyflakes gate"
    cat /tmp/pyflakes.out
fi
elapsed

# ── Phase 5: Grobid (CRF-only image) ──────────────────────────────
section "Phase 5 — Grobid"
# Clean stale container from a prior run.
sudo docker rm -f corpus-grobid 2>/dev/null || true
# CRF-only image is ~7 GB (vs ~32 GB for grobid/grobid:0.8.1) and
# produces the same REST surface for the demo's 11 PDFs.
sudo docker pull lfoppiano/grobid:0.8.1 > /tmp/grobid_pull.out 2>&1
# -XX:-UseContainerSupport: the JVM bundled with lfoppiano/grobid:0.8.1
# is old enough that its cgroup v2 detector hits a known NPE
# (CgroupV2Subsystem.getInstance → "anyController is null") on modern
# Ubuntu AMIs that default to cgroup v2. Disabling container-aware
# sizing skips the broken codepath; -Xmx/-Xms are already set
# explicitly so we don't actually need the auto-detection.
sudo docker run -d --name corpus-grobid -p 8070:8070 \
    -e JAVA_OPTS='-XX:-UseContainerSupport -Xmx4g -Xms1g' \
    --restart=unless-stopped \
    lfoppiano/grobid:0.8.1 > /dev/null

printf "  waiting for Grobid /api/isalive"
grobid_alive=0
for i in $(seq 1 120); do
    if [ "$(curl -fsS http://localhost:8070/api/isalive 2>/dev/null)" = "true" ]; then
        printf " — alive after %ds\n" "$((i * 2))"
        grobid_alive=1
        break
    fi
    sleep 2; printf '.'
done
if [ "$grobid_alive" -eq 1 ]; then
    note_pass "Grobid reachable at http://localhost:8070"
else
    note_fail "Grobid did not become ready within 240s"
    sudo docker logs --tail 50 corpus-grobid || true
    exit 1
fi
elapsed

# ── Phase 6: corpus check + corpus run ────────────────────────────
section "Phase 6 — corpus check + corpus run"
cd "$REPO_ROOT/demo"
rm -rf output  # clean slate so we exercise the full run, not --resume

if corpus check; then
    note_pass "corpus check exits 0"
else
    note_fail "corpus check failed (non-zero exit)"
fi

echo "  running 'corpus -v run --no-vision'; output → $LOG"
if corpus -v run --no-vision > "$LOG" 2>&1; then
    note_pass "corpus run completed (exit 0)"
else
    note_fail "corpus run failed; tail of $LOG follows"
    tail -60 "$LOG"
    exit 1
fi
elapsed

# ── Phase 7: programmatic verification ────────────────────────────
section "Phase 7 — verify success criteria"

# (a) corpus status: 11/11 across every stage row + no failures or flags.
status_out=$(corpus -v status --report 2>&1)
if echo "$status_out" | grep -q "Failures: none recorded" \
   && echo "$status_out" | grep -q "Quality flags: none recorded"; then
    note_pass "corpus status: no failures, no quality flags"
else
    note_fail "corpus status: failures or quality flags present"
    echo "$status_out" | grep -E "Failures:|Quality flags:|recorded" || true
fi
n_complete=$(echo "$status_out" | grep -cE "11 / 11" || true)
if [ "$n_complete" -ge 11 ]; then
    note_pass "corpus status: $n_complete stage rows at 11/11 (≥ 11)"
else
    note_fail "corpus status: only $n_complete stage rows at 11/11 (expected ≥ 11)"
fi

# (b) bundle_manifest.json shape.
manifest="$REPO_ROOT/demo/output/_serve/bundle_manifest.json"
if [ ! -f "$manifest" ]; then
    note_fail "bundle_manifest.json not written"
else
    paper_count=$(jq -r '.paper_count' "$manifest")
    chunk_count=$(jq -r '.chunk_count' "$manifest")
    figure_count=$(jq -r '.figure_count' "$manifest")
    bundle_version=$(jq -r '.bundle_version' "$manifest")
    [ "$paper_count" = "11" ] \
        && note_pass "manifest.paper_count = 11" \
        || note_fail "manifest.paper_count = $paper_count (expected 11)"
    [ "$chunk_count" -ge 800 ] \
        && note_pass "manifest.chunk_count = $chunk_count (≥ 800 — expect ~938 on macOS)" \
        || note_fail "manifest.chunk_count = $chunk_count (expected ≥ 800)"
    [ "$figure_count" -ge 100 ] \
        && note_pass "manifest.figure_count = $figure_count (≥ 100 — expect ~154 on macOS)" \
        || note_fail "manifest.figure_count = $figure_count (expected ≥ 100)"
    note_pass "manifest.bundle_version = $bundle_version"
fi

# (c) instructions.md shepherded (#4 — cli.py copies from corpuscle root).
[ -f "$REPO_ROOT/demo/output/instructions.md" ] \
    && note_pass "instructions.md present in output_dir root" \
    || note_fail "instructions.md missing from output_dir root"
[ -f "$REPO_ROOT/demo/output/_serve/instructions.md" ] \
    && note_pass "instructions.md present in served bundle" \
    || note_fail "instructions.md missing from served bundle"

# (d) #70 absolute-path audit clean.
if grep -qE "Path scrub: rewrote [0-9]+ files; audit clean\." "$LOG"; then
    note_pass "served-bundle absolute-path audit clean (#70)"
else
    note_fail "served-bundle audit did not log 'audit clean' — check $LOG"
fi

# (e) Known regressions we must NOT see in the run log.
# Each pattern is the "smoke ran cleanly" guarantee for a prior fix.
neg_patterns=(
    "Table 'document_chunks' already exists"     # #71
    "Could not apply staged bib overrides"       # 3393a48 (db_path NameError)
    "Expected top-level file instructions.md missing"  # a04aff6 (shepherd)
    "Dynamo is not supported on Python"          # macOS-Rosetta marker
    "Bad address.*g\\+\\+"                       # #72 (inductor JIT)
    "name 'db_path' is not defined"              # belt-and-suspenders for 3393a48
    "NameError: name '"                          # any new NameError class regression
)
for pat in "${neg_patterns[@]}"; do
    if grep -qE "$pat" "$LOG"; then
        note_fail "regression: run log contains \"$pat\""
    else
        note_pass "run log clean of \"$pat\""
    fi
done

# (f) MCP / SSE round-trip via the canonical smoke tool.
cd "$REPO_ROOT"
if python tools/smoke_test_sse.py demo/output/_serve --port 18080 \
       > /tmp/sse_smoke.out 2>&1; then
    note_pass "tools/smoke_test_sse.py — all four layers passed"
else
    note_fail "tools/smoke_test_sse.py failed; tail follows"
    tail -40 /tmp/sse_smoke.out
fi

elapsed

# ── Summary ────────────────────────────────────────────────────────
section "Summary"
printf "  passed: %s%d%s\n" "$GREEN" "$PASS_COUNT" "$RESET"
printf "  failed: %s%d%s\n" "$RED" "$FAIL_COUNT" "$RESET"

{
    echo "passed=$PASS_COUNT"
    echo "failed=$FAIL_COUNT"
    if [ "$FAIL_COUNT" -gt 0 ]; then
        printf 'fail_detail: %s\n' "${FAIL_DETAILS[@]}"
    fi
    echo "log=$LOG"
} > "$SUMMARY"

echo
if [ "$FAIL_COUNT" -eq 0 ]; then
    printf '%s═══ EC2 SMOKE PASSED ═══%s\n' "$GREEN" "$RESET"
    echo "Bundle: $REPO_ROOT/demo/output/_serve/"
    echo "Run log: $LOG"
    echo "Summary: $SUMMARY"
    exit 0
else
    printf '%s═══ EC2 SMOKE FAILED (%d criteria) ═══%s\n' \
        "$RED" "$FAIL_COUNT" "$RESET"
    echo "Failures:"
    for d in "${FAIL_DETAILS[@]}"; do echo "  - $d"; done
    echo "Full run log: $LOG"
    exit 1
fi

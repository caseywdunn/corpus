#!/usr/bin/env bash
# clean_install_walkthrough.sh
#
# End-to-end UX walkthrough: bootstrap a fresh conda env, install the
# corpus package, build a small corpuscle from a directory of PDFs,
# serve it over MCP, and exercise every operator-facing verb at least
# once (corpus init / check / run / status / serve / bib / completion,
# plus --dry-run / --help / --cite / --version).
#
# This is a *reference* — adapt the paths near the top to your own
# host and corpuscle. Authored while sanity-checking the v0.3 CLI
# surface; produced most of the UX fixes in the surrounding commit
# range. Re-run after CLI-affecting changes (new subcommands, log
# format tweaks, transport flags) to confirm the operator path stays
# clean.
#
# Run interactively, copy-paste section by section. Several blocks are
# meant to be inspected, not piped. Don't `bash` end-to-end.
#
# The default paths assume the worked example used to author this
# script (a 12-paper siphonophore corpus). Edit these three vars to
# point at your own scratch dir, the corpus repo, and a directory
# that carries your .bib / lexicon.yaml / instructions.md.
#
# Layout this script produces under EXERCISE_DIR:
#   ./pdfs/         10 short siphonophore PDFs (<20 pages each), the
#                   primary input for `corpus run`
#   ./pdfs_extra/   2 additional PDFs reserved for the second run
#                   (cp into ./pdfs/ at §8 to exercise resume)

set -euo pipefail   # only honored if you actually `bash` a section.

EXERCISE_DIR="$HOME/runs/corpus_exercise"
CORPUS_REPO="$HOME/repos/corpus"
SIPH_REPO="$HOME/repos/siphonophores"


# =============================================================================
# 0. Populate ./pdfs/ and ./pdfs_extra/ from the siphonophores library.
#    Idempotent — re-running re-copies the same 12 files, so the exercise
#    stays reproducible. Page counts in the comments are from
#    `pdfinfo <pdf> | grep ^Pages`; they all satisfy "<20 pages each".
# =============================================================================

mkdir -p "$EXERCISE_DIR/pdfs" "$EXERCISE_DIR/pdfs_extra"
cd "$EXERCISE_DIR"

# 10 PDFs for the first run.
LIB="$SIPH_REPO/library"
cp "$LIB/D_E/Daniel_Daniel1963b.pdf"      pdfs/   #  3 pages
cp "$LIB/F/Fewkes1889f.pdf"               pdfs/   # 10 pages
cp "$LIB/G/Grossmannetal2013.pdf"         pdfs/   #  8 pages
cp "$LIB/H_J/Iinuma_etal2020.pdf"         pdfs/   # 13 pages
cp "$LIB/K/Kramp1943.pdf"                 pdfs/   # 19 pages
cp "$LIB/K/Krsinic1993.pdf"               pdfs/   #  4 pages
cp "$LIB/L/Lo_etal2013.pdf"               pdfs/   # 15 pages
cp "$LIB/M/Mackie1962.pdf"                pdfs/   #  3 pages
cp "$LIB/M/Mackie2002.pdf"                pdfs/   #  5 pages
cp "$LIB/U_Z/Williams_Conway1981.pdf"     pdfs/   # 12 pages

# 2 PDFs for the second run (added at §8 to exercise implicit resume).
cp "$LIB/A/Alvarino1980b.pdf"             pdfs_extra/   #  6 pages
cp "$LIB/P/Pugh_Gasca2009.pdf"            pdfs_extra/   #  8 pages

ls pdfs/        | wc -l       # should print 10
ls pdfs_extra/  | wc -l       # should print 2

# Verify each PDF is short enough (sanity check, no failures expected).
for f in pdfs/*.pdf pdfs_extra/*.pdf; do
  pages=$(pdfinfo "$f" 2>/dev/null | awk '/^Pages:/ {print $2}')
  printf '%3s  %s\n' "$pages" "$(basename "$f")"
done


# =============================================================================
# 1. Clean conda env (named corpus_test so it does not collide with the
#    dev `corpus` env). Mirrors the README "Installation" block.
# =============================================================================

# If a stale corpus_test already exists, blow it away (safe — nothing else
# uses the name).
conda env remove -n corpus_test -y 2>/dev/null || true

# environment.yaml defaults `name: corpus`. Override with -n so the env is
# named corpus_test instead.
conda env create -n corpus_test -f "$CORPUS_REPO/environment.yaml"

# Activate. `conda activate` only works after sourcing the conda hook in a
# non-interactive shell — drop the eval line if your shell already has it.
eval "$(conda shell.bash hook)"
conda activate corpus_test

# Editable install of the corpus package (puts the `corpus` binary on PATH).
pip install -e "$CORPUS_REPO"

# Sanity-check the entry point + version.
corpus --version
corpus --help


# =============================================================================
# 2. Tesseract language packs (needed for OCR on older scans).
# =============================================================================

bash "$CORPUS_REPO/tools/install_tessdata.sh"


# =============================================================================
# 3. Grobid (PDF metadata + reference parsing). Runs as a docker service;
#    `corpus run` does NOT auto-launch it.
# =============================================================================

# docker-compose.yml lives in the corpus repo, so launch from there.
( cd "$CORPUS_REPO" && docker compose up -d grobid )

# Probe — should print "true" once warm (first start can take ~30s while
# DeLFT models load).
until curl -sf http://localhost:8070/api/isalive | grep -q true; do
  echo "waiting for grobid..."; sleep 3
done


# =============================================================================
# 4. Scaffold the corpuscle. We are already in $EXERCISE_DIR with ./pdfs/
#    populated; `corpus init` writes config.yaml from the bundled template.
# =============================================================================

cd "$EXERCISE_DIR"

corpus init                   # writes ./config.yaml

# The template defaults are already close to what we want:
#   input_pdfs: ./pdfs
#   taxonomy.source: worms / root_id: 1267
# We need to point bib + lexicon at the curated files in ~/repos/siphonophores.
# Easiest: overwrite the two commented lines with absolute paths.
python - <<PY
from pathlib import Path
p = Path("config.yaml")
text = p.read_text()
text = text.replace(
    "# bib: ./references.bib",
    f"bib: {Path.home()}/repos/siphonophores/siphonophores.bib",
)
text = text.replace(
    "# lexicon: ./lexicon.yaml",
    f"lexicon: {Path.home()}/repos/siphonophores/lexicon.yaml",
)
p.write_text(text)
print(text)
PY

# Drop instructions.md alongside config.yaml so `corpus serve` picks it up
# automatically (server reads <output_dir>/instructions.md by default;
# also accepts `corpus serve --instructions <path>`).
cp "$SIPH_REPO/instructions.md" .


# =============================================================================
# 5. Pre-flight checks.
# =============================================================================

corpus check                  # config schema + grobid + GPU + disk
corpus run --dry-run          # prints the plan; no artifacts written


# =============================================================================
# 6. First full pipeline run (10 PDFs).
#
# CPU-only laptop: expect ~10–30 minutes for 10 short PDFs the first time
# (taxonomy.sqlite from WoRMS dominates the wall clock on a cold cache).
# Re-runs are idempotent — only changed inputs do work.
# =============================================================================

corpus run

# Postmortem rollup.
corpus status --report
corpus status --json | head -40
corpus status --sort-by quality_flag_count --tail 5 || true


# =============================================================================
# 7. Local MCP server smoke-test.
#
# Two ways to validate:
#   (a) `corpus serve --check` — non-binding pre-flight (bundle present,
#       DBs loadable, port free when paired with --transport sse). Exits 0
#       on green.
#   (b) Bring it up over SSE on a local port and probe with curl. The
#       default stdio transport is meant for an MCP client launching the
#       process, so SSE is what you actually test from a shell.
#
# Port choice: we deliberately avoid 8080 because code-server (browser
# VS Code) ships its own service there on EC2 dev hosts; corpus serve
# would silently fail to bind and curl would hit code-server instead,
# producing a confusing `{"error":"Unauthorized"}` response from the
# wrong server. Pick something out of the way; sanity-check with:
#   ss -tlnp | grep -E ":<port> "
# (no output → free)
MCP_PORT=8090

# (a) Pre-flight only. --transport sse triggers the port-free check;
#     without it the bundle/DB checks still run but the port isn't probed.
corpus serve -- --check --transport sse --port "$MCP_PORT"

# (b) SSE smoke test. Bearer-auth is required for SSE; generate a throwaway
#     token, start the server in the background, hit /sse, then stop it.
echo "test-token-$(date +%s)" > .mcp-token
chmod 600 .mcp-token

corpus serve -- \
    --transport sse --host 127.0.0.1 --port "$MCP_PORT" \
    --auth-token-file .mcp-token \
    > serve.log 2>&1 &
SERVE_PID=$!

# Wait for the port. Probe /healthz (non-streaming, returns "ok"
# immediately) — NOT /sse, which holds the connection open as a server-
# sent-event stream and would hang `curl -sf` indefinitely without a
# --max-time. Cap the loop at ~40s for cold-start index loading.
for i in {1..40}; do
  curl -sf -o /dev/null --max-time 1 \
    "http://127.0.0.1:$MCP_PORT/healthz" && break
  sleep 1
done

# Hit the SSE endpoint, dump the first event, then move on. A healthy
# server answers with `event: endpoint` followed by a session URL
# (sample output: `event: endpoint\ndata: /messages/?session_id=...`).
# A 401 means you're talking to a different server on this port (see
# note above about code-server on 8080).
curl -N --max-time 3 \
    -H "Authorization: Bearer $(cat .mcp-token)" \
    "http://127.0.0.1:$MCP_PORT/sse" || true

# Sanity-check auth rejection (both should 401).
echo "--- wrong token ---"
curl -s -o /dev/null -w "HTTP %{http_code}\n" --max-time 3 \
    -H "Authorization: Bearer wrongtoken" \
    "http://127.0.0.1:$MCP_PORT/sse"
echo "--- no auth ---"
curl -s -o /dev/null -w "HTTP %{http_code}\n" --max-time 3 \
    "http://127.0.0.1:$MCP_PORT/sse"

# Bring it down.
kill "$SERVE_PID" 2>/dev/null || true
wait "$SERVE_PID" 2>/dev/null || true
tail -20 serve.log


# =============================================================================
# 8. Add two more PDFs and re-run. Implicit resume only re-processes what
#    changed; the 10 already-processed papers should skip every stage.
# =============================================================================

cp pdfs_extra/*.pdf pdfs/
ls pdfs/ | wc -l       # should be 12

corpus run --dry-run

corpus run             # idempotent; new PDFs only

corpus status --report


# =============================================================================
# 9. Re-test the local server against the updated bundle.
# =============================================================================

corpus serve -- --check --transport sse --port "$MCP_PORT"

corpus serve -- \
    --transport sse --host 127.0.0.1 --port "$MCP_PORT" \
    --auth-token-file .mcp-token \
    > serve2.log 2>&1 &
SERVE_PID=$!
# Wait via /healthz (same reasoning as §7).
for i in {1..40}; do
  curl -sf -o /dev/null --max-time 1 \
    "http://127.0.0.1:$MCP_PORT/healthz" && break
  sleep 1
done
curl -N --max-time 3 \
    -H "Authorization: Bearer $(cat .mcp-token)" \
    "http://127.0.0.1:$MCP_PORT/sse" || true
kill "$SERVE_PID" 2>/dev/null || true
wait "$SERVE_PID" 2>/dev/null || true


# =============================================================================
# 10. Tour the rest of the CLI surface: help, cite, dry-run variants,
#     bib round-trip, completion, taxonomy module.
# =============================================================================

# Top-level + per-verb help.
corpus --help
corpus run --help
corpus serve --help
corpus status --help
corpus bib --help
corpus check --help
corpus completion --help

# Citation (CITATION.cff is the source of truth).
corpus --cite               # plain text
corpus --cite=bibtex        # BibTeX

# Version (includes a git short-SHA when sourced from a dev clone).
corpus --version

# Verbosity knobs.
corpus -v run --dry-run     # INFO logging
corpus -vv run --dry-run    # DEBUG logging
corpus -q run --dry-run     # WARNING+ only

# Plan-only variants of run (no artifacts written).
corpus run --dry-run                       # full plan
corpus run --dry-run --no-vision           # skip vision pass
corpus run --dry-run --no-bundle           # skip served-bundle distillation
corpus run --dry-run --no-prune            # audit orphans only, don't delete
corpus run --dry-run --skip-checks         # bypass grobid/input_pdfs preflight
corpus run --dry-run --force-rebuild       # would rebuild every cross-paper DB

# Cross-paper-DB selective rebuilds (each is a real run, so only paste if you
# want them to actually execute).
# corpus run --force-rebuild-taxonomy
# corpus run --force-rebuild-biblio
# corpus run --force-rebuild-taxon-mentions

# Bib round-trip (export → edit → import). Source of truth for user edits.
corpus bib export -o my_corpus.bib
ls -la my_corpus.bib
# $EDITOR my_corpus.bib                    # fix titles, years, authors
# corpus bib import my_corpus.bib          # re-apply edits

# Status filters and triage.
corpus status --report
corpus status --json | python -m json.tool | head -60
corpus status --sort-by quality_flag_count --tail 10
corpus status --sort-by stage_failure_count --tail 10
corpus status --propose-skips --min-flags 2
corpus status --skipped
corpus status --list-hashes --filter-gate empty_text || true

# Shell completion.
corpus completion bash > /tmp/corpus.bash
corpus completion zsh  > /tmp/corpus.zsh
corpus completion fish > /tmp/corpus.fish
head -20 /tmp/corpus.bash

# Direct module entry points (corpus subcommands shell out to these).
python -m pipeline.taxonomy_ingest --help
python -m mcpsrv.bundle --help
python -m bib.export --help
python -m bib.importer --help


# =============================================================================
# 11. Tear down.
# =============================================================================

# Stop grobid (frees ~6 GB of RAM).
( cd "$CORPUS_REPO" && docker compose stop grobid )

# Drop the env entirely if you're done with this exercise.
# conda deactivate
# conda env remove -n corpus_test -y

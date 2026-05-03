#!/usr/bin/env bash
# Wrapper that picks the right corpus-env Python for this host and execs
# mcp_server.py. Invoked from .mcp.json so the same committed config works
# on a local Mac and on Bouchet.

set -e

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"

case "$OSTYPE" in
    darwin*)
        # Apple Silicon needs a native arm64 env (Intel anaconda cannot
        # provide a torch new enough for current docling/transformers).
        CORPUS_PY="${CORPUS_PY:-$HOME/miniforge3/envs/corpus/bin/python}"
        # conda-forge numpy and torch each ship their own libomp.dylib;
        # set this so both can coexist in one process.
        export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"
        ;;
    linux*)  CORPUS_PY="${CORPUS_PY:-/home/cwd7/.conda/envs/corpus/bin/python}" ;;
    *)       echo "run_mcp_server.sh: unknown OSTYPE=$OSTYPE" >&2; exit 1 ;;
esac

if [[ ! -x "$CORPUS_PY" ]]; then
    echo "run_mcp_server.sh: $CORPUS_PY not found or not executable" >&2
    echo "Override with CORPUS_PY=/path/to/python" >&2
    exit 1
fi

OUTPUT_DIR="${CORPUS_OUTPUT_DIR:-$REPO_DIR/output}"

cd "$REPO_DIR"
exec "$CORPUS_PY" mcp_server.py "$OUTPUT_DIR"

#!/usr/bin/env bash
# install_ocr_extras.sh — install ocrmypdf's optional image-compression
# helpers into the active environment.
#
# Why this exists: without `pngquant`, pipeline/scan.py degrades
# ocrmypdf from `--optimize 2` to `--optimize 1`, and on scanned
# historical material that is not a minor size difference. Measured on
# Beklemishev 1969 (a 45-page Russian scan) with the v1.0 re-OCR path:
#
#     source PDF                     1.1 MB
#     --optimize 1 (no pngquant)      90 MB
#     --optimize 2 (with pngquant)    35 MB     <- 61% smaller
#
# Across a whole library that is the difference between a corpuscle that
# fits a quota and one that doesn't. `corpus run` re-OCRs every scanned
# paper since v1.0, so far more documents take this path than used to.
#
# Usage:
#     conda activate corpus
#     bash tools/install_ocr_extras.sh
#
# Idempotent: anything already on PATH is left alone.
#
# Not in environment.yaml because conda env files have no platform
# conditionals and `pngquant` is absent from conda-forge osx-arm64,
# which is a supported target — adding it there would break `conda env
# create` on Apple Silicon. `jbig2enc` is not on conda-forge for any
# platform, so it stays a system install everywhere.

set -euo pipefail

have() { command -v "$1" >/dev/null 2>&1; }

os="$(uname -s)"
arch="$(uname -m)"

echo "Installing ocrmypdf image-compression helpers"
echo "Platform: $os/$arch"
echo

# --- pngquant -------------------------------------------------------
# Colour-PNG quantization. Required by ocrmypdf for --optimize 2+.
if have pngquant; then
    printf '  %-10s already on PATH (%s)\n' "pngquant" "$(command -v pngquant)"
elif [[ -n "${CONDA_PREFIX:-}" ]] && [[ "$os/$arch" != "Darwin/arm64" ]]; then
    printf '  %-10s installing from conda-forge...\n' "pngquant"
    conda install -y -q -c conda-forge pngquant
    have pngquant && printf '  %-10s ok (%s)\n' "pngquant" "$(command -v pngquant)"
else
    printf '  %-10s NOT INSTALLED\n' "pngquant"
    if [[ "$os" == "Darwin" ]]; then
        echo "       conda-forge has no osx-arm64 build. Install with:"
        echo "           brew install pngquant"
    else
        echo "       Activate the conda env first, or install at system level:"
        echo "           sudo apt-get install pngquant"
    fi
fi

# --- jbig2enc -------------------------------------------------------
# Bitonal (black-and-white) image compression at the highest optimization
# level. Not packaged on conda-forge for any platform.
if have jbig2; then
    printf '  %-10s already on PATH (%s)\n' "jbig2enc" "$(command -v jbig2)"
else
    printf '  %-10s NOT INSTALLED (no conda-forge package on any platform)\n' "jbig2enc"
    case "$os" in
        Darwin) echo "       brew install jbig2enc" ;;
        Linux)  echo "       sudo apt-get install jbig2enc   # or build from source" ;;
    esac
    echo "       On HPC without root, check first:  module avail jbig2enc"
    echo "       Optional — ocrmypdf runs without it, at a lower compression"
    echo "       level on bitonal scans."
fi

echo
echo "Current state:"
for tool in pngquant jbig2; do
    if have "$tool"; then
        printf '  %-10s yes\n' "$tool"
    else
        printf '  %-10s no\n' "$tool"
    fi
done
echo
if have pngquant; then
    echo "ocrmypdf will run at the configured ocr.optimize_level (default 2)."
else
    echo "ocrmypdf will be degraded to --optimize 1; expect substantially"
    echo "larger processed.pdf files on scanned material."
fi

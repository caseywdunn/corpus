#!/usr/bin/env bash
# install_tessdata.sh — download Tesseract LSTM language packs into the
# active conda env's tessdata/ directory.
#
# Why this exists: conda-forge's `tesseract` package ships only the
# English LSTM model.  None of the `tesseract-data-<code>` packages
# referenced in older copies of environment.yaml actually exist on
# conda-forge — they're a different upstream's naming.  The portable
# install path is to drop `<code>.traineddata` files from the official
# tesseract-ocr/tessdata_best repo into $CONDA_PREFIX/share/tessdata,
# which is exactly what tesseract searches at runtime.
#
# Usage:
#     conda activate corpus
#     bash tools/install_tessdata.sh                # default language set
#     bash tools/install_tessdata.sh ara hin tha    # explicit list
#
# Pip-only / non-conda installs:
#     TESSDATA_DIR=/usr/local/share/tessdata bash tools/install_tessdata.sh
#
# Idempotent: existing files are skipped.  All downloads come from
# tessdata_best (highest-accuracy LSTM models) for parity with what
# the pipeline was tuned against on historical scans.

set -euo pipefail

# Default fallback set + opt-in classical packs, matching README.md
# §"Language support".  deu_latf (19th-c. German Fraktur) is bundled
# in the default set — the v0.1 docs treated it as a separate manual
# install, but the curl pattern is identical, so combining is cleaner.
DEFAULT_LANGS=(
    deu fra rus lat spa por
    chi_sim chi_tra jpn ell kor
    grc          # Ancient Greek — opt-in, but small and harmless
    deu_latf     # 19th-c. German Fraktur (Goldfuss, Pagenstecher, Brandt, Dönitz)
)

BASE_URL="https://github.com/tesseract-ocr/tessdata_best/raw/main"

if [[ $# -gt 0 ]]; then
    LANGS=("$@")
else
    LANGS=("${DEFAULT_LANGS[@]}")
fi

# Resolve target dir.  Conda env first; explicit override wins for
# users who installed tesseract via brew/apt instead.
if [[ -n "${TESSDATA_DIR:-}" ]]; then
    TESSDATA="$TESSDATA_DIR"
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
    TESSDATA="$CONDA_PREFIX/share/tessdata"
else
    echo "ERROR: neither TESSDATA_DIR nor CONDA_PREFIX is set." >&2
    echo "       Activate the corpus conda env first:" >&2
    echo "           conda activate corpus" >&2
    echo "       Or set TESSDATA_DIR to your tesseract data directory." >&2
    exit 1
fi

if [[ ! -d "$TESSDATA" ]]; then
    echo "ERROR: tessdata directory not found at $TESSDATA." >&2
    echo "       Is the 'tesseract' package installed in this env?" >&2
    exit 1
fi

echo "Installing Tesseract language packs into $TESSDATA"
echo "Source: tessdata_best (high-accuracy LSTM)"
echo

for lang in "${LANGS[@]}"; do
    target="$TESSDATA/$lang.traineddata"
    if [[ -f "$target" ]]; then
        printf '  %-10s already present, skipping\n' "$lang"
        continue
    fi
    printf '  %-10s downloading... ' "$lang"
    # --retry covers transient GitHub raw-content hiccups; --fail makes
    # 4xx/5xx exit non-zero so set -e catches them.
    if curl -fsSL --retry 3 --retry-delay 2 \
            -o "$target.tmp" "$BASE_URL/$lang.traineddata"; then
        mv "$target.tmp" "$target"
        size=$(wc -c < "$target" | tr -d ' ')
        printf 'ok (%d bytes)\n' "$size"
    else
        rm -f "$target.tmp"
        echo "FAILED" >&2
        echo "       URL: $BASE_URL/$lang.traineddata" >&2
        echo "       (typo in language code? check the ISO->Tesseract map" >&2
        echo "        in pipeline/scan.py)" >&2
        exit 1
    fi
done

echo
echo "Installed packs in $TESSDATA:"
ls -1 "$TESSDATA"/*.traineddata 2>/dev/null \
    | sed 's|.*/||;s|\.traineddata$||' \
    | sort \
    | column || true

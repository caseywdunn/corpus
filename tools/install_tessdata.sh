#!/usr/bin/env bash
# install_tessdata.sh — download Tesseract LSTM language packs into the
# active conda env's tessdata/ directory.
#
# Why this exists: the pipeline is tuned against tessdata_best (the
# highest-accuracy LSTM models), and none of the `tesseract-data-<code>`
# packages referenced in older copies of environment.yaml actually exist
# on conda-forge — they're a different upstream's naming.  The portable
# install path is to drop `<code>.traineddata` files from the official
# tesseract-ocr/tessdata_best repo into $CONDA_PREFIX/share/tessdata,
# which is exactly what tesseract searches at runtime.
#
# conda-forge's `tesseract` (>= 5.5) bundles ~158 packs of its own, but
# they are **tessdata_fast** builds — byte-identical to
# tesseract-ocr/tessdata_fast, the smallest and least accurate of the
# three upstream variants (e.g. rus is 3.9 MB fast vs 15.3 MB best).
# A plain "file exists → skip" therefore silently left the fast models
# in place on every fresh conda env and downloaded only the one code
# conda-forge does not ship (deu_latf), which is the opposite of what
# the banner below claims.  We track what *this script* installed in a
# marker file and re-fetch anything else, so the packs the pipeline
# actually uses are always tessdata_best.
#
# Usage:
#     conda activate corpus
#     bash tools/install_tessdata.sh                # default language set
#     bash tools/install_tessdata.sh ara hin tha    # explicit list
#     bash tools/install_tessdata.sh --force        # re-download everything
#
# Pip-only / non-conda installs:
#     TESSDATA_DIR=/usr/local/share/tessdata bash tools/install_tessdata.sh
#
# Idempotent: packs this script has already installed are skipped.

set -euo pipefail

# Default fallback set + opt-in classical packs, matching README.md
# §"Language support".  deu_latf (19th-c. German Fraktur) is bundled
# in the default set — the v0.1 docs treated it as a separate manual
# install, but the curl pattern is identical, so combining is cleaner.
DEFAULT_LANGS=(
    eng          # always appended as the last OCR fallback, so it matters
                 # most; conda-forge ships only a tessdata_fast build of it
    deu fra rus lat spa por
    chi_sim chi_tra jpn ell kor
    grc          # Ancient Greek — opt-in, but small and harmless
    deu_latf     # 19th-c. German Fraktur (Goldfuss, Pagenstecher, Brandt, Dönitz)
)

BASE_URL="https://github.com/tesseract-ocr/tessdata_best/raw/main"

FORCE=0
REQUESTED=()
for arg in "$@"; do
    case "$arg" in
        --force|-f) FORCE=1 ;;
        -h|--help)
            sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        -*)
            echo "ERROR: unknown option '$arg'." >&2
            echo "       Usage: bash tools/install_tessdata.sh [--force] [lang ...]" >&2
            exit 1
            ;;
        *) REQUESTED+=("$arg") ;;
    esac
done

if [[ ${#REQUESTED[@]} -gt 0 ]]; then
    LANGS=("${REQUESTED[@]}")
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

# Records the codes this script has fetched from tessdata_best.  A pack
# present on disk but absent from the marker came from somewhere else
# (almost always the conda-forge tesseract package's tessdata_fast
# bundle) and is replaced.
MARKER="$TESSDATA/.corpus_tessdata_best"
touch "$MARKER"

echo "Installing Tesseract language packs into $TESSDATA"
echo "Source: tessdata_best (high-accuracy LSTM)"
echo

for lang in "${LANGS[@]}"; do
    target="$TESSDATA/$lang.traineddata"
    if [[ -f "$target" ]] && [[ $FORCE -eq 0 ]] \
            && grep -qxF "$lang" "$MARKER"; then
        printf '  %-10s already present, skipping\n' "$lang"
        continue
    fi
    if [[ -f "$target" ]] && ! grep -qxF "$lang" "$MARKER"; then
        printf '  %-10s replacing non-best copy... ' "$lang"
    else
        printf '  %-10s downloading... ' "$lang"
    fi
    # --retry covers transient GitHub raw-content hiccups; --fail makes
    # 4xx/5xx exit non-zero so set -e catches them.
    if curl -fsSL --retry 3 --retry-delay 2 \
            -o "$target.tmp" "$BASE_URL/$lang.traineddata"; then
        mv -f "$target.tmp" "$target"
        grep -qxF "$lang" "$MARKER" || echo "$lang" >> "$MARKER"
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

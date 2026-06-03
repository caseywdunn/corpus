# Single source of truth for Bouchet paths.
# Sourced by every batch_*.sh script in this repo, and referenced by
# BOUCHET.md for one-time setup commands.
#
# Override any of these by exporting before sbatch'ing, e.g.:
#     OUTPUT_DIR=$BOUCHET_PROJECT/output_sample sbatch batch_pass3b.sh

# Lab project storage (PI Casey Dunn). All persistent state lives under
# this path on Bouchet — never use $HOME for corpus data, since $HOME has
# a 125 GB quota that the model cache + LanceDB index will exceed.
BOUCHET_PROJECT="${BOUCHET_PROJECT:-/nfs/roberts/project/pi_cwd7/cwd7}"

REPO_DIR="${REPO_DIR:-$BOUCHET_PROJECT/corpus}"

# Corpuscle config.yaml — the source of truth for the new `corpus run`
# flow (#138). input_pdfs / output_dir / bib / lexicon / taxonomy / grobid
# now live INSIDE this file, not as CLI flags. The phase scripts
# (batch_process_corpus, batch_pass3b, batch_embed, batch_finalize) drive
# `corpus -c "$CORPUS_CONFIG" run --only <phase>`. See BOUCHET.md for how
# to author it. The dynamic Grobid node URL is still injected at submit
# time via $GROBID_URL, which `corpus run` honors as an override.
# Corpuscles live under $BOUCHET_PROJECT/corpuscles/<name>_YYYYMMDD/.
# Update this line when starting a new build; the date is the build date.
# Sample runs follow the same convention (siphonophore_sample_YYYYMMDD).
CORPUS_CONFIG="${CORPUS_CONFIG:-$BOUCHET_PROJECT/corpuscles/siphonophore_20260603/config.yaml}"

# Legacy path vars — still consumed by the not-yet-ported helper scripts
# (batch_grobid.sh, batch_biblio.sh) and by batch_pipeline.sh's banner.
# Keep them in sync with the matching keys in $CORPUS_CONFIG.
INPUT_DIR="${INPUT_DIR:-$BOUCHET_PROJECT/siphonophores/library}"
OUTPUT_DIR="${OUTPUT_DIR:-$BOUCHET_PROJECT/corpuscles/siphonophore_20260603}"
BIB_FILE="${BIB_FILE:-$BOUCHET_PROJECT/siphonophores/siphonophores.bib}"
LEXICON="${LEXICON:-$BOUCHET_PROJECT/siphonophores/lexicon.yaml}"
CACHE_DIR="${CACHE_DIR:-$BOUCHET_PROJECT/cache}"

# HuggingFace cache (BGE-M3 + Qwen2.5-VL weights). Pre-download once on
# a compute node — see BOUCHET.md.
export HF_HOME="${HF_HOME:-$CACHE_DIR/huggingface}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME/hub}"

# Native libs for the conda env (#20). Bouchet compute nodes don't
# always have $CONDA_PREFIX/lib on the runtime loader path; without
# it docling/torch silently fall back to whatever the system has and
# produce mis-structured output without crashing. Prepend the env's
# libs so they win — matches the resolution order pip / conda used
# at install time. Caller must `conda activate corpus` before sourcing.
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
fi

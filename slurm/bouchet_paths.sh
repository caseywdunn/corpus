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
INPUT_DIR="${INPUT_DIR:-$BOUCHET_PROJECT/siphonophores}"
OUTPUT_DIR="${OUTPUT_DIR:-$BOUCHET_PROJECT/output}"
CACHE_DIR="${CACHE_DIR:-$BOUCHET_PROJECT/cache}"

# HuggingFace cache (BGE-M3 + Qwen2.5-VL weights). Pre-download once on
# a compute node — see BOUCHET.md.
export HF_HOME="${HF_HOME:-$CACHE_DIR/huggingface}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME/hub}"

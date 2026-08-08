# Installation notes

The one-command conda install in the [README](README.md) covers most setups. The additions below are needed only for specific OCR or deployment scenarios.

## How corpus is distributed

**git clone + conda is the canonical channel, and the only supported one.** Pin a release with `git clone --branch v1.0.0 https://github.com/caseywdunn/corpus.git`.

**There is no PyPI package, and this is deliberate rather than pending.** corpus depends on native tools pip cannot install — tesseract, ghostscript, pngquant, pandoc, and a Grobid server — so `pip install corpus` would produce a package that imports cleanly and then fails on the first scanned PDF. That is precisely the failure class the 1.0 cycle existed to remove, and shipping a channel that reintroduces it would be a step backwards. The [pip-only fallback](#pip-only-fallback) below exists for hosts that already have the native toolchain; it is a fallback, not a distribution channel.

conda-forge *could* carry the native toolchain, but a feedstock needs every dependency on conda-forge, and the pip-only block (`docling`, `lancedb`, `mcp`, `ocrmypdf`, …) is not. A container image is the better future channel — Docker is already a prerequisite for Grobid — and is tracked in [PLAN.md](dev_docs/PLAN.md) rather than promised here.

**For citation, use the Zenodo DOI**, which is what a PyPI presence would have been standing in for: [10.5281/zenodo.19964909](https://doi.org/10.5281/zenodo.19964909). Each release gets its own DOI via the GitHub integration, and [CITATION.cff](CITATION.cff) carries the current form.

## Supported platforms

| Target | Status |
| --- | --- |
| **linux-x86_64** (HPC clusters, generic CPU/GPU servers) | Supported |
| **macOS arm64** (Apple Silicon) | Supported — see [Apple Silicon: arm64-native conda required](#apple-silicon-arm64-native-conda-required) below |
| macOS x86_64 (Intel Mac, or Rosetta on Apple Silicon) | Not supported. Apple dropped Intel-mac PyTorch wheels after 2.2; `docling` and `transformers ≥ 5.0` both require torch ≥ 2.4. The chain is structurally broken. |
| linux-aarch64 | Not currently supported. Most wheels exist, but `pymupdf` (among others) isn't on conda-forge linux-aarch64, so `environment.yaml` would need pip-side workarounds. Add as a third target if a real Graviton use case appears. |

## Apple Silicon: arm64-native conda required

On an arm64 Mac the corpus environment must be created with an arm64-native conda. An x86_64 conda on Apple Silicon runs under Rosetta and traps the install in the unsupported macOS Intel matrix above. Symptoms: `transformers` silently disables PyTorch, then `docling` model loads abort with `Dynamo is not supported on Python 3.12+`.

Anaconda, Miniconda, and Miniforge all ship arm64 builds for Apple Silicon, and their downloads pages typically pick the right one for an Apple Silicon visitor. The hazard is an older install that predates Apple Silicon support, picking the wrong installer from a multi-button downloads page, or running conda from inside a Rosetta'd terminal. Confirm via `conda info | grep platform` (see the pre-creation gate below). Distribution download pages:

- **Anaconda** — <https://www.anaconda.com/download> (Apple Silicon Graphical Installer).
- **Miniconda** — <https://docs.conda.io/projects/miniconda/en/latest/> (`Miniconda3-latest-MacOSX-arm64.sh`).
- **Miniforge** — <https://github.com/conda-forge/miniforge> (arm64 by default since 2021).

### Pre-creation gate

Check whichever conda you intend to use:

```bash
conda info | grep platform
# must print: platform : osx-arm64   (NOT osx-64)
```

If it prints `osx-arm64`, you're set — proceed to `conda env create -f environment.yaml` and skip to the post-creation gate.

If it prints `osx-64`, your conda is running under Rosetta. **Recommended fix:** install one of the arm64 distributions above alongside the existing conda. Use its explicit binary path so the env lands under the right distribution; the existing conda is left untouched.

```bash
# Example with miniforge — substitute any of the three arm64 distributions
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
~/miniforge3/bin/conda env create -f environment.yaml
```

**Alternative for users who specifically don't want a second conda distribution:** `CONDA_SUBDIR=osx-arm64` can force the existing x86_64 conda to pull arm64 packages.

```bash
CONDA_SUBDIR=osx-arm64 conda env create -f environment.yaml
conda activate corpus && conda config --env --set subdir osx-arm64
```

The second line locks `subdir: osx-arm64` into the env's `.condarc` so subsequent `conda update` calls don't revert to x86_64. Less battle-tested than the primary path — the conda binary itself still runs under Rosetta, and a missing subdir-lock will silently re-pull x86_64 packages on update.

### Post-creation gate

After `conda env create`, verify the env's Python is arm64-native. `corpus check` enforces this automatically (it hard-fails on Rosetta'd Python on macOS), but the manual equivalent is:

```bash
conda activate corpus
python -c "import platform; print(platform.machine())"
# expect: arm64   (NOT x86_64)
```

`platform.machine()` reports the architecture of the running process, not the silicon — so it reads `arm64` only when the env's Python is truly native, not when an x86_64 binary is being launched under Rosetta on Apple Silicon. Running it from the env's Python is the meaningful check; running it from `/usr/bin/python3` or another unrelated Python tells you nothing about the corpus env.

`pip install -e .` then puts the `corpus` binary inside `<conda-prefix>/envs/corpus/bin/`. Activate the env or call the binary by absolute path — Claude Desktop / VS Code MCP configs that previously pointed at `/opt/anaconda3/envs/corpus/bin/corpus` need to be updated to wherever the arm64 conda put the env.

## Model downloads and where they land

Four models are fetched from HuggingFace the first time you run the
pipeline, not at install time: docling's page-layout model (~165 MB),
docling's TableFormer (~340 MB), docling's `HybridChunker` tokenizer
(`sentence-transformers/all-MiniLM-L6-v2`, ~730 KB), and
[BGE-M3](https://huggingface.co/BAAI/bge-m3) (~4.3 GB) for embeddings.
A fifth, the ~16 GB Qwen2.5-VL local vision model, is fetched only for
`--figure-panels vision-local`.

The chunker tokenizer is tiny and easy to overlook, but it is not
optional: without it `HybridChunker` fails and chunking silently falls
back to a naive character window — a two-page paper came out as 1 chunk
instead of 16.

```bash
corpus prefetch                    # get them now, not mid-run
corpus prefetch --include-vision   # ...plus the local VLM
```

Prefetching retries with backoff, because HuggingFace throttles anonymous
traffic with HTTP 429 — a shared institutional NAT looks like abuse from
the other side ([#140](https://github.com/caseywdunn/corpus/issues/140)).
`corpus check` reports what is cached **without touching the network**,
so it is safe on an isolated host.

They land in `~/.cache/huggingface/` unless you set `HF_HOME` (or the
narrower `HF_HUB_CACHE`) *before* downloading — do that when your home
directory is small or isn't shared across the machines that will run the
pipeline.

### Offline hosts

Prepare where there is internet, run where there isn't:

```bash
export HF_HOME=$PROJECT/huggingface   # shared storage, both places
corpus prefetch                       # where there is outbound access

export HF_HUB_OFFLINE=1               # where there isn't
corpus run
```

`HF_HUB_OFFLINE=1` is the important half: without it a missing model
stalls instead of failing. It is also worth setting on a connected host,
where it pins the run to the snapshot you cached instead of whatever
upstream `main` points at today. The taxonomy needs the same treatment —
`source: worms` walks a REST API, so pre-build it once and switch the
corpuscle to `source: dwca`
([#139](https://github.com/caseywdunn/corpus/issues/139)).

## Clusters with a small or non-writable home directory

On clusters where home space is unusable and work lives in a shared
project filesystem, conda, pip, and HuggingFace all need redirecting
*before* anything is installed, or they quietly fill `$HOME` and fail
([#153](https://github.com/caseywdunn/corpus/issues/153)).

```bash
export CONDA_PKGS_DIRS=$PROJECT/conda/pkgs
export CONDA_ENVS_DIRS=$PROJECT/conda/envs
export PIP_CACHE_DIR=$SCRATCH/pip-cache
export HF_HOME=$PROJECT/huggingface
export TMPDIR=$SCRATCH/tmp

conda env create -f environment.yaml && conda activate corpus
pip install -e .
bash tools/install_tessdata.sh
bash tools/install_ocr_extras.sh
corpus prefetch
```

Some sites additionally need `$HOME` itself redirected for the install,
because not every tool honors the variables above — prefer the specific
variables, and reach for `HOME=$PROJECT/fakehome conda env create …` only
if something still writes to your real home.

Grobid runs as a batch job rather than a Docker service on a cluster, and
`grobid.url` must then name the allocated compute node rather than
`localhost`. [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) is the worked
example, including the Singularity invocation and SLURM scripts.

## Higher OCR compression: pngquant + jbig2enc

`ocrmypdf` has two optional native helpers. The runtime auto-degrades when they're missing (`pipeline/scan.py` drops `--optimize` from 2 → 1 when `pngquant` isn't on PATH), so this is a size-of-output decision — but on scanned material it is a large one, and since v1.0 re-OCRs every scan rather than trusting its text layer, most of a historical corpus now takes that path.

Measured on Beklemishev 1969, a 45-page Russian scan:

| | output |
|---|---|
| source PDF | 1.1 MB |
| `--optimize 1` (no pngquant) | 90 MB |
| `--optimize 2` (with pngquant) | **35 MB** |

Across a library that is the difference between fitting a quota and not. The helper script installs what your platform has:

```bash
conda activate corpus
bash tools/install_ocr_extras.sh
```

It is idempotent and prints what it could not install. The manual equivalents follow.

`pngquant` — needed for color-PNG quantization at `--optimize 2+`. On conda-forge for linux-64 and osx-64 but **not osx-arm64**, which is why it can't live in `environment.yaml` (conda env files have no platform conditionals, so listing it would break `conda env create` on Apple Silicon). `tools/install_ocr_extras.sh` installs it from conda-forge where available; otherwise:

```bash
# macOS
brew install pngquant

# Debian/Ubuntu
sudo apt-get install pngquant
```

`jbig2enc` — B/W image compression at the highest optimization level. Not packaged on conda-forge for any platform, so it is a system install everywhere:

```bash
# macOS
brew install jbig2enc

# Debian/Ubuntu
sudo apt-get install jbig2enc
```

**On an HPC cluster without root**, `bash tools/install_ocr_extras.sh` still gets you `pngquant`: it installs from conda-forge into the active env, which your account owns. Run it on a login node — `conda install` needs outbound network, which compute nodes frequently lack.

`jbig2enc` is the one you may have to skip. It has no conda package on any platform, and `apt` / `brew` both need root. It is also the smaller win — bitonal images only, while `pngquant` covers the colour and greyscale scans that dominate historical literature. Check `module avail jbig2enc`; if it isn't there, skipping it is a reasonable default rather than a compromise, and `ocrmypdf` falls back gracefully.

## OCR language packs

The pipeline is tuned against [`tesseract-ocr/tessdata_best`](https://github.com/tesseract-ocr/tessdata_best), the highest-accuracy of the three upstream LSTM variants. Those `<code>.traineddata` files have to be dropped into `$CONDA_PREFIX/share/tessdata/` by hand: there is no `tesseract-data-<code>` package on conda-forge; older versions of `environment.yaml` referenced packages by that name and silently failed to install on a fresh env (issue [#52](https://github.com/caseywdunn/corpus/issues/52)).

The conda-forge `tesseract` package *does* bundle ~158 packs of its own, so a fresh env looks fully equipped — but they are [`tessdata_fast`](https://github.com/tesseract-ocr/tessdata_fast) builds, byte-identical to that upstream and the least accurate variant (`rus` 3.9 MB vs 15.3 MB for best; `deu` 1.5 MB vs 8.6 MB). Only `deu_latf`, the 19th-century German Fraktur pack, is genuinely absent. Running the installer is what upgrades the rest.

[`tools/install_tessdata.sh`](tools/install_tessdata.sh) automates the download for the default fallback set in `config.yaml` (`eng`, `deu`, `fra`, `rus`, `lat`, `spa`, `por`, `chi_sim`, `chi_tra`, `jpn`, `ell`, `kor`, `grc`, plus `deu_latf` for 19th-c. German):

```bash
conda activate corpus
bash tools/install_tessdata.sh
bash tools/install_tessdata.sh --force   # re-download even the ones it installed
```

The script is idempotent: it records the codes it fetched in `$CONDA_PREFIX/share/tessdata/.corpus_tessdata_best` and skips those on a re-run. A pack that is on disk but *not* in that marker came from the conda-forge bundle and gets replaced — that is the whole point, so don't be surprised when a fresh env reports `replacing non-best copy` for every language.

To add a language outside the default set, pass its [ISO-to-Tesseract code](pipeline/scan.py) explicitly:

```bash
bash tools/install_tessdata.sh ara hin tha
```

For non-conda installs (`pip install -r requirements.txt` + a system-installed tesseract), point the script at your tessdata directory directly:

```bash
TESSDATA_DIR=/usr/local/share/tessdata bash tools/install_tessdata.sh
# Debian/Ubuntu typically: /usr/share/tesseract-ocr/<version>/tessdata
```

## Grobid image choice, and hosts without Docker

The default is the lightweight `lfoppiano/grobid:0.8.1` (~1 GB, CRF-only)
— the same image used on HPC, so local and cluster match. The full DeLFT
image `grobid/grobid:0.8.1` (~32 GB) parses references better but requires
AVX, so it needs an AVX-capable Linux x86_64 host; it crash-loops under
Rosetta on Apple Silicon. Opt in via the comments in
[`docker-compose.yml`](docker-compose.yml).

Without Docker, [Singularity/Apptainer](https://docs.sylabs.io/) pulls the
same image: `singularity build grobid.sif docker://lfoppiano/grobid:0.8.1`.
On a cluster it runs as a job, so `grobid.url` must name the allocated node
— see [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) for the full recipe.

## Pip-only fallback

If you can't use conda, you'll need to install the system tools yourself (`brew install ghostscript tesseract pngquant jbig2enc pandoc` on macOS, or `apt-get install` the equivalents on Debian/Ubuntu) and then:

```bash
pip install -e ".[dev]"   # development clone (the [dev] extra adds
                          # pytest, pyflakes, and ipykernel — they are
                          # not runtime dependencies, see #162)
pip install -e .          # runtime only, e.g. a server-only deploy host

# or for a deploy host pinning a release:
pip install git+https://github.com/caseywdunn/corpus.git@v1.1.0
```

Python **3.12** specifically: `requires-python` is `>=3.12,<3.13`, because
that is what `environment.yaml` pins and what CI tests, and the pinned ML
stack (torch 2.12, transformers 5.8.1, docling 2.94.0) has been verified
against nothing else. `python3 -m venv` on a current distro may well give
you 3.13, so create the venv with an explicit interpreter:

```bash
python3.12 -m venv .venv && . .venv/bin/activate
```

`requirements.txt` is retained for the AWS deploy path (which builds a stock `python3.12 -m venv` and pre-existed `pyproject.toml`); both manifests must stay in sync per [CONTRIBUTING.md](CONTRIBUTING.md#dependencies--two-files-on-purpose).

## Updating an existing environment

After `environment.yaml` changes:

```bash
conda env update -f environment.yaml --prune
```

If a previous attempt left a partial `corpus` environment behind, remove it first with `conda env remove -n corpus`.

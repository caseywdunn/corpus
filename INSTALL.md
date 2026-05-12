# Installation notes

The one-command conda install in the [README](README.md) covers most setups. The additions below are needed only for specific OCR or deployment scenarios.

## Supported platforms

| Target | Status |
| --- | --- |
| **linux-x86_64** (HPC clusters, generic CPU/GPU servers, Bouchet) | Supported |
| **macOS arm64** (Apple Silicon) | Supported — see [Apple Silicon: use miniforge, not Intel anaconda](#apple-silicon-use-miniforge-not-intel-anaconda) below |
| macOS x86_64 (Intel Mac, or Rosetta on Apple Silicon) | **Not supported.** Apple dropped Intel-mac PyTorch wheels after 2.2; `docling` and `transformers ≥ 5.0` both require torch ≥ 2.4. The chain is structurally broken. |
| linux-aarch64 | Not currently supported. Most wheels exist, but `pymupdf` (among others) isn't on conda-forge linux-aarch64, so `environment.yaml` would need pip-side workarounds. Add as a third target if a real Graviton use case appears. |

## Apple Silicon: use miniforge, not Intel anaconda

On an arm64 Mac the corpus environment must be created with an arm64-native conda. The default download from anaconda.com is the **Intel x86_64** build of Anaconda3, which keeps running under Rosetta and traps the install in the unsupported macOS Intel matrix above. Symptoms when this happens: `transformers` silently disables PyTorch, then `docling` model loads abort with `Dynamo is not supported on Python 3.12+`.

[Miniforge](https://github.com/conda-forge/miniforge) ships native arm64 builds. Install it once:

```bash
# https://github.com/conda-forge/miniforge#install
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
# then use the miniforge conda (not /opt/anaconda3) for `conda env create`
~/miniforge3/bin/conda env create -f environment.yaml
```

Verify the env is actually arm64 after creation:

```bash
~/miniforge3/envs/corpus/bin/python -c "import platform; print(platform.machine())"
# expect: arm64   (NOT x86_64)
```

`pip install -e .` then puts the `corpus` binary inside `~/miniforge3/envs/corpus/bin/`. Activate the env or call the binary by absolute path — Claude Desktop / VS Code MCP configs that previously pointed at `/opt/anaconda3/envs/corpus/bin/corpus` need to be updated.

## Higher OCR compression: pngquant + jbig2enc

`ocrmypdf` has two optional native helpers, neither of which is available on conda-forge for every platform we target. The runtime auto-degrades when they're missing (`pipeline/scan.py` drops `--optimize` from 2 → 1 when `pngquant` isn't on PATH), so installing them is purely a size-of-output decision.

`pngquant` — needed for color-PNG quantization at `--optimize 2+`. Missing from conda-forge **osx-arm64**, so we install it at the system level on all platforms for consistency:

```bash
# macOS
brew install pngquant

# Debian/Ubuntu
sudo apt-get install pngquant
```

`jbig2enc` — needed for B/W image compression at the highest optimization level:

```bash
# macOS
brew install jbig2enc

# Debian/Ubuntu
sudo apt-get install jbig2enc
```

On Bouchet, `module avail pngquant jbig2enc` will tell you whether modules are available; if not, skip — `ocrmypdf` falls back gracefully.

## OCR language packs

The conda-forge `tesseract` package ships only the English LSTM model. Every other language pack — including 19th-century German Fraktur (`deu_latf`) — has to be dropped into `$CONDA_PREFIX/share/tessdata/` as a `<code>.traineddata` file from the [`tesseract-ocr/tessdata_best`](https://github.com/tesseract-ocr/tessdata_best) repo. There is no `tesseract-data-<code>` package on conda-forge; older versions of `environment.yaml` referenced packages by that name and silently failed to install on a fresh env (issue [#52](https://github.com/caseywdunn/corpus/issues/52)).

[`tools/install_tessdata.sh`](tools/install_tessdata.sh) automates the download for the default fallback set in `config.yaml` (`eng`, `deu`, `fra`, `rus`, `lat`, `spa`, `por`, `chi_sim`, `chi_tra`, `jpn`, `ell`, `kor`, `grc`, plus `deu_latf` for 19th-c. German):

```bash
conda activate corpus
bash tools/install_tessdata.sh
```

The script is idempotent — re-running skips packs that already exist.

To add a language outside the default set, pass its [ISO-to-Tesseract code](pipeline/scan.py) explicitly:

```bash
bash tools/install_tessdata.sh ara hin tha
```

For non-conda installs (`pip install -r requirements.txt` + a system-installed tesseract), point the script at your tessdata directory directly:

```bash
TESSDATA_DIR=/usr/local/share/tessdata bash tools/install_tessdata.sh
# Debian/Ubuntu typically: /usr/share/tesseract-ocr/<version>/tessdata
```

## Grobid on Bouchet

Grobid runs as a Docker service locally. On an HPC cluster without Docker, [Singularity](https://docs.sylabs.io/) can pull the same image:

```bash
singularity build grobid.sif docker://lfoppiano/grobid:0.8.1
singularity run --bind $HOME grobid.sif &
corpus run    # config.yaml in cwd points at <input> + sets grobid.url
```

## Pip-only fallback

If you can't use conda, you'll need to install the system tools yourself (`brew install ghostscript tesseract pngquant jbig2enc pandoc` on macOS, or `apt-get install` the equivalents on Debian/Ubuntu) and then:

```bash
pip install -e .          # development clone
# or for a deploy host pinning a release:
pip install git+https://github.com/caseywdunn/corpus.git@v0.3.0
```

`requirements.txt` is retained for the AWS deploy path (which builds a stock `python3.12 -m venv` and pre-existed `pyproject.toml`); both manifests must stay in sync per [CONTRIBUTING.md](CONTRIBUTING.md#dependencies--two-files-on-purpose).

## Updating an existing environment

After `environment.yaml` changes:

```bash
conda env update -f environment.yaml --prune
```

If a previous attempt left a partial `corpus` environment behind, remove it first with `conda env remove -n corpus`.

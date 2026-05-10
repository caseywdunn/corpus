# Installation notes

The one-command conda install in the [README](README.md) covers most setups. The additions below are needed only for specific OCR or deployment scenarios.

## Higher OCR compression: jbig2enc

`pngquant` ships with the conda environment. `jbig2enc` is the other optional `ocrmypdf` helper and is not packaged on conda-forge for all platforms — install it at the system level for further size savings:

```bash
# macOS
brew install jbig2enc

# Debian/Ubuntu
sudo apt-get install jbig2enc
```

On Bouchet, `module avail jbig2enc` will tell you whether a module is available; if not, skip — `ocrmypdf` falls back gracefully.

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
python process_corpus.py <input> <output> --grobid-url http://localhost:8070
```

## Pip-only fallback

If you can't use conda, you'll need to install the system tools yourself (`brew install ghostscript tesseract pngquant jbig2enc` on macOS, or `apt-get install` the equivalents on Debian/Ubuntu) and then:

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

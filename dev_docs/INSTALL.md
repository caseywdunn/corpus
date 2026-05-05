# Installation notes

The one-command conda install in the [README](../README.md) covers most setups. The additions below are needed only for specific OCR or deployment scenarios.

## Higher OCR compression: jbig2enc

`pngquant` ships with the conda environment. `jbig2enc` is the other optional `ocrmypdf` helper and is not packaged on conda-forge for all platforms — install it at the system level for further size savings:

```bash
# macOS
brew install jbig2enc

# Debian/Ubuntu
sudo apt-get install jbig2enc
```

On Bouchet, `module avail jbig2enc` will tell you whether a module is available; if not, skip — `ocrmypdf` falls back gracefully.

## Additional OCR language packs

The base `tesseract` install includes only English. For a multilingual / historical corpus, add the language data you need:

```bash
# Modern European languages + Latin
conda install -c conda-forge tesseract-data-deu tesseract-data-fra tesseract-data-rus tesseract-data-lat
```

**19th-century German Fraktur** (e.g., Haeckel, Schneider) is not on conda-forge. Download `deu_latf.traineddata` from the Tesseract project and drop it into the env's `tessdata` directory:

```bash
TESSDATA="$CONDA_PREFIX/share/tessdata"
curl -L -o "$TESSDATA/deu_latf.traineddata" \
  https://github.com/tesseract-ocr/tessdata_best/raw/main/deu_latf.traineddata
```

## Grobid on Bouchet

Grobid runs as a Docker service locally. On an HPC cluster without Docker, [Singularity](https://docs.sylabs.io/) can pull the same image:

```bash
singularity build grobid.sif docker://lfoppiano/grobid:0.8.1
singularity run --bind $HOME grobid.sif &
python process_corpus.py <input> <output> --grobid-url http://localhost:8070
```

## Pip-only fallback

If you can't use conda, `requirements.txt` still works but you'll need to install the system tools yourself (`brew install ghostscript tesseract pngquant jbig2enc` on macOS, or `apt-get install` the equivalents on Debian/Ubuntu):

```bash
pip install -r requirements.txt
```

## Updating an existing environment

After `environment.yaml` changes:

```bash
conda env update -f environment.yaml --prune
```

If a previous attempt left a partial `corpus` environment behind, remove it first with `conda env remove -n corpus`.

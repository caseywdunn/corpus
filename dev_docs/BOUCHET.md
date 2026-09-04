# Running the corpus pipeline on Bouchet (YCRC)

Operational notes for running the full siphonophore corpus on Bouchet (~2,000 PDFs).
The library grows as papers are added, so treat every corpus-size figure here as
approximate. For the current count, run `corpus check` (it reports `input_pdfs`
with a total) or count directly — note the library is nested in letter
subdirectories, so this needs `find`, not `ls`:

```bash
find "$BOUCHET_PROJECT/siphonophores/library" -iname '*.pdf' | wc -l
```

See [OVERVIEW.md](OVERVIEW.md) for pipeline architecture and
`dunnlab-hpc` skill for Dunn-lab-wide YCRC conventions.

## One-time setup

### 1. Repo + corpus on project storage

All paths in this guide and the `batch_*.sh` scripts derive from a single env var, `$BOUCHET_PROJECT`, defined in [slurm/bouchet_paths.sh](../slurm/bouchet_paths.sh). The default is `/nfs/roberts/project/pi_cwd7/cwd7` (PI Casey Dunn's lab project storage). To run from a different project root, edit that one line — every batch script sources it.

```bash
export BOUCHET_PROJECT=/nfs/roberts/project/pi_cwd7/cwd7
cd "$BOUCHET_PROJECT"
git clone git@github.com:caseywdunn/corpus.git
git clone git@github.com:dunnlab/siphonophores.git

# git-lfs is not on Bouchet's default PATH — load the module first.
# If `module load git-lfs` fails, search for the right name with
# `module spider lfs`, or install into the conda env as a fallback:
#   conda install -n corpus -c conda-forge git-lfs
module load git-lfs
cd siphonophores
git lfs install --local                   # wires LFS hooks into .git/config
git lfs pull                              # ~2,000 PDFs, several GB
```

Resulting layout (all under `$BOUCHET_PROJECT`):

```
corpus/                          ← this repo (the `corpus` CLI)
siphonophores/                   ← input PDF tree + siphonophores.bib + lexicon.yaml
corpuscles/                      ← all siphonophore corpuscle builds
  current -> siphonophore_YYYYMMDD  ← symlink naming the active build
  siphonophore_YYYYMMDD/         ← one directory per production build (date = build date)
    config.yaml                  ← authored in step 3 (the source of truth)
    documents/<HASH>/…           ← per-paper artifacts (created by extract)
    *.sqlite, vector_db/         ← cross-paper DBs + LanceDB (created by embed/post)
    _serve/                      ← distilled served bundle (created by bundle)
  siphonophore_gold_YYYYMMDD/    ← smoke-test builds over the 35 transcribed
                                   documents (same structure; see below)
cache/huggingface/               ← model cache (see below)
cache/grobid.sif                 ← Grobid Singularity image (step 6)
cache/grobid_{tmp,logs}/<jobid>/ ← per-job Grobid scratch; each submit sweeps
                                   dirs left by earlier jobs (step 6)
```

When starting a new build, create `corpuscles/siphonophore_YYYYMMDD/`, scaffold a
`config.yaml` there (step 3), and repoint the `current` symlink at it:

```bash
cd "$BOUCHET_PROJECT/corpuscles"
ln -sfn siphonophore_YYYYMMDD "$BOUCHET_PROJECT/corpuscles/current"
ls -l "$BOUCHET_PROJECT/corpuscles/current"   # confirm the active build
```

The link *target* is a bare name, resolved relative to the link's own directory,
so this command is correct from any cwd — the `cd` is for the reader, not the shell.

**No tracked file names a corpuscle**, so switching builds is never a commit.
The symlink lives on the cluster next to the builds it points at, which also
makes "which corpuscle are jobs writing to?" answerable with `ls` rather than
`git log`.

As of #138 the SLURM phase scripts drive the **same `corpus run` CLI**
users run, one phase per job. All per-corpuscle inputs (PDF dir, BibTeX,
lexicon, taxonomy source, Grobid) live in the corpuscle's `config.yaml`
— **not** as CLI flags or env vars. The scripts reference it through
`$CORPUS_CONFIG`, resolved in
[slurm/bouchet_paths.sh](../slurm/bouchet_paths.sh) as: exported
`$CORPUS_CONFIG` → `corpuscles/current/config.yaml` → hard error. There is
deliberately no dated fallback; a missing selection fails loudly rather than
silently resuming a finished build. `corpus` honors `$CORPUS_CONFIG` natively
too (#61), and `bouchet_paths.sh` exports the resolved value so every job in a
submitted chain builds the same corpuscle even if `current` moves mid-run.

To build a corpuscle other than `current` for one run, export the var — no
symlink change, no edit:

```bash
cd "$BOUCHET_PROJECT/corpus"
CORPUS_CONFIG=$BOUCHET_PROJECT/corpuscles/siphonophore_gold_YYYYMMDD/config.yaml \
    bash slurm/batch_pipeline.sh
```

### 2. Conda environment

```bash
module load miniconda
conda env create -f "$BOUCHET_PROJECT/corpus/environment.yaml"
conda activate corpus
pip install -e "$BOUCHET_PROJECT/corpus"
bash "$BOUCHET_PROJECT/corpus/tools/install_tessdata.sh"
bash "$BOUCHET_PROJECT/corpus/tools/install_ocr_extras.sh"
```

**Don't skip `install_tessdata.sh`.** The conda-forge `tesseract` package bundles
~158 packs of its own, so a fresh env *looks* fully equipped — but they are
low-accuracy [`tessdata_fast`](https://github.com/tesseract-ocr/tessdata_fast)
builds (`rus` 3.9 MB against 15.3 MB for `tessdata_best`), and no `tessdata_best`
pack is installable from conda-forge (issue #52). The script downloads the default
`ocr.ocr_languages_default` set (deu, fra, rus, lat, spa, por, chi_sim, chi_tra,
jpn, ell, kor, grc, plus 19th-c. German Fraktur `deu_latf`) into
`$CONDA_PREFIX/share/tessdata/`. It is idempotent, so re-running is safe. Skip it
and every non-English paper silently OCRs as English — including the Fraktur
scans the smoke test below is meant to exercise. To add a language outside the
default set, pass its code: `bash tools/install_tessdata.sh ara hin tha`. See
[INSTALL.md](../INSTALL.md#ocr-language-packs).

**`install_ocr_extras.sh` matters more here than on a laptop.** It installs
`pngquant` from conda-forge into the env — no root needed, but run it on a login
node since `conda install` wants outbound network. Without it ocrmypdf drops from
`--optimize 2` to `--optimize 1`, and since v1.0 re-OCRs every scan rather than
trusting its text layer, that is a large disk difference across a full library: a
45-page Russian scan came out 90 MB instead of 35 MB. `jbig2enc` has no conda
package anywhere and needs root, so check `module avail jbig2enc` and skip it
otherwise — it only affects bitonal images.

Confirm before moving on:

```bash
ls "$CONDA_PREFIX/share/tessdata/deu_latf.traineddata"   # must exist
```

The step 7 preflight re-verifies this properly — `corpus check` reports every
configured pack, e.g. `all 13 configured packs installed`.

**GPU torch (issue #21) — no extra step needed.** The conda env is complete as
created: `environment.yaml` pins `torch==2.12.0`, PyPI ships that as a CUDA 13
build, and Bouchet's GPU driver (580.159.04 / CUDA 13.0) runs it natively.

**Do not `pip install` torch by hand here.** `torch==2.12.0` does not exist on
PyTorch's `cu128` index (it stops at 2.11.0), so a cu128 re-install silently
downgrades torch and breaks the #98 known-good pin. If a future driver rollback
ever makes a CUDA 12 build necessary, 2.12.0 is available from the `cu126` and
`cu130` indexes — use one of those and keep the pin.

Verify on a **GPU node** (torch can't see a device from the login node):

```bash
srun -p gpu_devel --gpus=h200:1 -t 5 -c 2 --mem=8G \
    python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
# → 2.12.0+cu130 True
```

Use `gpu_devel` for interactive GPU checks. `srun`/`salloc` against `gpu`,
`gpu_h200`, `gpu_b200` or `gpu_rtx6000` fails with `Invalid qos specification`:
those partitions' `AllowQos` is `normal,nothrottle` only. That is partition
policy, not a broken allocation, and it doesn't affect `sbatch` — which is how
every pipeline phase is submitted. `day`, `devel`, `gpu_devel` and `scavenge` do
allow interactive jobs.

The Stage 1 batch script's preflight (`slurm/batch_process_corpus.sh`) aborts
loudly if torch / docling can't import, so a botched env fails fast instead of
producing mis-structured output.

*(Verified on the cluster 2026-08-01.)*

### 3. Author the corpuscle config.yaml

This is the source of truth for the build — `corpus run` reads every
per-corpuscle input from it. Scaffold it once, then edit:

```bash
conda activate corpus
# Replace YYYYMMDD with today's date, e.g. 20260603
mkdir -p "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD"
cd "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD"
corpus init                              # drops a commented config.yaml here
```

Repoint the `current` symlink at the new directory (see step 1) so the batch
scripts pick it up:

```bash
cd "$BOUCHET_PROJECT/corpuscles"
ln -sfn siphonophore_YYYYMMDD "$BOUCHET_PROJECT/corpuscles/current"
ls -l "$BOUCHET_PROJECT/corpuscles/current"   # confirm the active build
```

(As in step 1: the bare-name target makes this correct from any cwd; the `cd` is
for the reader.)

Edit `config.yaml` so it points at the siphonophores repo. Paths resolve
**relative to the config file's directory** (here, `corpuscles/siphonophore_YYYYMMDD/`),
so the `../../siphonophores/...` paths below work regardless of
`$BOUCHET_PROJECT` (YAML has no env-var expansion — keep them relative or
write absolute literals):

```yaml
# $BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/config.yaml
input_pdfs: ../../siphonophores/library       # the PDF tree
output_dir: .                                 # artifacts land here in the
                                              #   corpuscle dir;
                                              #   keeps resume pointed at any
                                              #   existing documents/<HASH>/ tree
bib: ../../siphonophores/siphonophores.bib
lexicon: ../../siphonophores/lexicon.yaml

taxonomy:                                     # Siphonophorae, WoRMS AphiaID 1371
  source: dwca
  path: ../../siphonophores/taxonomy.dwca.zip # pre-exported; see step 4
  root_id: 1371

figures:
  panel_detection: vision-local               # Pass 3b uses local Qwen2.5-VL

grobid:
  url: http://localhost:8070                  # overridden at submit time by
                                              #   $GROBID_URL (see step 6 / #138)
```

Don't run `corpus check` yet — steps 4–6 supply the taxonomy DB, the model cache,
and Grobid, so a check run here would fail on prerequisites you haven't built.
Step 7 is the preflight gate, and by then it should come back green.

Once `current` points at this corpuscle, the phase scripts find it with no
extra flags. To build a *different* corpuscle for one run, export
`CORPUS_CONFIG=/path/to/other/config.yaml` before submitting — it takes
precedence over the symlink.

### 4. Pre-build the taxonomy

**Required, not an optimization — do it before submitting any batch jobs.**
The array tasks run `corpus run --only extract`, and that path does *not*
build `taxonomy.sqlite` the way a full `corpus run` does: the orchestrator
hard-errors ~60 s in, and `afterok` then takes Pass 3b, Embed and Finalize
down with it. Two consecutive siphonophore builds lost their whole chain this
way (#251). `slurm/batch_pipeline.sh` now pre-builds it for you before the
first `sbatch`, so following the runbook is enough — but if you submit
`batch_process_corpus.sh` by hand, this step is yours.

The second reason is the original one: you do not want one extract task per
paper each walking the WoRMS REST API. Build `taxonomy.sqlite` once, up
front — it reads the source and path straight from `config.yaml`, and
no-ops when the file already exists.

**If this corpuscle was built with 1.2.1**, run it once on 1.2.2 or later:
that release's automatic pre-build doubled the `names` table on every launch
(#262), and the next ingest deduplicates it in place and logs how many rows
it removed. Lookups were correct throughout; the table was just growing.

```bash
corpus taxonomy ingest      # reads taxonomy.{source,path,root_id} from config.yaml
```

or spell the source out explicitly:

```bash
conda activate corpus
corpus -c "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/config.yaml" taxonomy ingest --source worms --root-id 1371
# → writes corpuscles/siphonophore_YYYYMMDD/taxonomy.sqlite (~700 KB, takes ~1–2 min)
```

**If `taxonomy.dwca.zip` is already committed to the siphonophores repo** (the
normal case after the first build), skip the WoRMS step and ingest directly
from the zip — works on any node, takes under a second:

```bash
conda activate corpus
corpus -c "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/config.yaml" \
    taxonomy ingest --source dwca \
    --input "$BOUCHET_PROJECT/siphonophores/taxonomy.dwca.zip" \
    --root-id 1371
# → writes corpuscles/siphonophore_YYYYMMDD/taxonomy.sqlite in ~1 s
```

Then export it as a portable DwC-A zip and commit it to the siphonophores repo.
This makes every future build ingest from the zip in seconds instead of
re-walking the WoRMS REST API, and keeps the taxonomy reproducible and
version-pinned rather than tracking a moving upstream:

```bash
corpus -c "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/config.yaml" taxonomy export \
    -o "$BOUCHET_PROJECT/siphonophores/taxonomy.dwca.zip"
```

Update `config.yaml` (the corpuscle dir is two levels below `siphonophores/`):

```yaml
taxonomy:
  source: dwca
  path: ../../siphonophores/taxonomy.dwca.zip  # local; no network, version-pinned
  root_id: 1371
```

Commit `taxonomy.dwca.zip` to the siphonophores repo so the taxonomy is
version-controlled alongside the PDFs and bib file. If WoRMS upstream changes
significantly, regenerate with `--source worms` and re-export.

### 5. Pre-download HuggingFace models

Do this once, before the first submission. Run it wherever is convenient — the
login node is fine (it pulls from the HuggingFace CDN at ~100 MB/s), and so are
batch nodes. Observe the usual courtesy: a multi-GB fetch on the shared login
node is fine, a long CPU-bound job is not.

Prefetching matters for two reasons: it keeps a full-corpus run from stopping
mid-stage on a HuggingFace 429, and it is what makes `HF_HUB_OFFLINE=1` usable
as a reproducibility lever (see below).

```bash
# Pre-download into the SAME cache the batch jobs read —
# $BOUCHET_PROJECT/cache/huggingface, the path slurm/bouchet_paths.sh
# exports as HF_HOME. Don't use ~/project/... here: if it doesn't
# resolve to $BOUCHET_PROJECT, the jobs won't find these weights and
# re-download them (the "Stale HF_HOME" pitfall below).
export HF_HOME="$BOUCHET_PROJECT/cache/huggingface"

# Everything the pipeline needs: docling's page-layout model and
# TableFormer, plus BGE-M3 for embeddings (~4.8 GB total).
corpus prefetch

# Add the Qwen2.5-VL-7B local vision backend (~15 GB) only if you will
# run Pass 3b with `--figure-panels vision-local`.
corpus prefetch --include-vision
```

Use `corpus prefetch` rather than hand-written `python -c` loader snippets: it
goes through the pipeline's own code paths, so it warms the full model set
(including docling's page-layout model and TableFormer, which the loader
one-liners missed) and picks up `HF_HUB_DISABLE_IMPLICIT_TOKEN=1` (#97) without
you exporting it. It also retries with backoff on HTTP 429.

`HF_HOME` is the only cache variable that matters. `TRANSFORMERS_CACHE` is dead
— `transformers` 5.x ignores it entirely — so don't set it.

`corpus check` reports what is already cached **without touching the
network**, so it is safe to run on a batch node to confirm the cache is
visible from there.

Batch jobs read the same `HF_HOME` (set by `slurm/bouchet_paths.sh` to
`$BOUCHET_PROJECT/cache/huggingface`), so they pick up these weights
without re-downloading.

**Is the pre-download idempotent?** *Partly.* Re-running `corpus prefetch`
is safe and a near-no-op when nothing changed upstream — it reports what is
already cached and skips it, and `huggingface_hub` HEAD-checks each file's
etag against the cache and fetches only what differs. But the loaders it
calls resolve the model's **`main` revision**, not a pinned commit, so the
download is
**not version-locked**: if `BAAI/bge-m3` or the Qwen repo publishes a new
revision upstream, the next run pulls it into a new cache snapshot. That
can change behavior — a changed bge-m3 silently makes fresh query vectors
inconsistent with a `vector_db/lancedb` index built from the old weights
(same dimension, different values), and a changed Qwen alters Pass 3b
output.

Check what's cached and which commit each model currently resolves to:

```bash
hf cache ls                                          # repos, REVISION (commit), refs, size, path
cat "$HF_HOME/hub/models--BAAI--bge-m3/refs/main"    # the commit `main` points at right now
```

For a reproducible build, set `HF_HUB_OFFLINE=1` in the batch environment
so a job can never reach out and pull a newer revision: it uses exactly the
snapshot cached here, or fails loudly. That is the lever that actually
holds, and it needs no code change. (Pinning a `revision="<sha>"` per
loader would be stricter still, but `corpus` does not currently thread a
revision through to the model loads — worth an issue if a build ever needs
to be reproducible against a moving upstream.)

### 6. Grobid as a Singularity service

Grobid is the extract phase's bibliographic + section-structure extractor. The
`docker-compose.yml` in this repo targets macOS dev; there is no Docker on Bouchet,
so run it via Singularity.

Build the image once (it already exists at `$BOUCHET_PROJECT/cache/grobid.sif`, so
skip this unless you're starting a fresh project root):

```bash
cd "$BOUCHET_PROJECT/cache"
singularity build grobid.sif docker://lfoppiano/grobid:0.8.1
```

**Start it with `slurm/batch_grobid.sh`, not by hand.** That script does three binds
that all matter — `$BOUCHET_PROJECT` (so Grobid can read the PDF tree), a writable
host dir onto `/opt/grobid/grobid-home/tmp` (without it Grobid answers **HTTP 500 to
every request**, because the Singularity rootfs is read-only and it writes temp files
per request), and a writable host dir onto `/opt/grobid/logs` (without it Grobid
**crashes on startup** opening `logs/grobid-service.log`). It creates both dirs
per-job so concurrent instances don't collide. A hand-rolled
`singularity run --bind "$BOUCHET_PROJECT"` does not even get that far: without
`--pwd /opt/grobid` it dies immediately with

```
[FATAL tini] exec ./grobid-service/bin/grobid-service failed: No such file or directory
```

because Singularity starts in your *host* cwd rather than the image's `WORKDIR`.
Add `--pwd` and it starts, then fails every request for want of the writable binds.

To bring one up and point this shell at it — this is what step 7 and a manual
`sbatch slurm/batch_process_corpus.sh` need:

```bash
cd "$BOUCHET_PROJECT/corpus"
GROBID_JOB=$(sbatch --parsable slurm/batch_grobid.sh)
until [ "$(squeue -j "$GROBID_JOB" -h -o %T)" = RUNNING ]; do sleep 5; done
export GROBID_URL="http://$(squeue -j "$GROBID_JOB" -h -o %N):8070"

# The job reaching RUNNING only means SLURM started the container; Grobid
# itself needs another ~30-60 s to load its models and bind :8070. The first
# few polls failing is expected, not an error — this is the same wait
# slurm/batch_pipeline.sh does for you.
for i in $(seq 1 60); do
    curl -fsS "$GROBID_URL/api/isalive" >/dev/null 2>&1 && break
    printf "\r  waiting for Grobid at %s (%ds)..." "$GROBID_URL" "$((i * 5))"
    sleep 5
done
echo
curl -fsS "$GROBID_URL/api/isalive"; echo   # → true
```

If that last `curl` errors instead of printing `true`, Grobid did not come up
within 5 minutes and the log is the place to look —
`logs/slurm-grobid-$GROBID_JOB.err` for model loading, and
`$BOUCHET_PROJECT/cache/grobid_logs/$GROBID_JOB/grobid-service.log` for the
service itself. A healthy startup ends with
`Started application@...{0.0.0.0:8070}`. Don't `scancel` a job that is merely
still loading Wapiti models.

Read that service log **while the job is alive** — the job deletes its contents
on exit. The now-empty `cache/grobid_logs/$GROBID_JOB/` directory usually
survives (NFS holds the still-open log as `.nfsXXXX` past the window SLURM
gives the cleanup trap), so leftover per-job directories under `cache/` are
expected rather than a symptom; the next `batch_grobid.sh` submit sweeps them.

`$GROBID_URL` overrides the config's `grobid.url` for this dynamically-allocated node
(#138), for both `corpus run` and `corpus check`. It is honored only when set, so
`config.yaml` stays authoritative on a standalone submit. Remember to `scancel
"$GROBID_JOB"` when you're done — it holds a `week` allocation for 48 h.

In practice you don't submit Grobid by hand — `slurm/batch_pipeline.sh`
(see [Production run](#production-run)) submits `slurm/batch_grobid.sh`, discovers the
node, waits for `/api/isalive`, exports `$GROBID_URL` into the extract job, and tears
Grobid down afterward. The manual path above exists for step 7 and for debugging —
and because the orchestrator always starts its own, **cancel a hand-started Grobid
before launching a production run** or it holds a `week` allocation for 48 h doing
nothing.

### 7. Preflight: `corpus check`

With steps 1–6 done, this is the gate. Run it from **inside the corpuscle** — config
resolution falls back to `./config.yaml`, so no `-c` and no long quoted path:

```bash
cd "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD"
corpus check
corpus run --dry-run
```

`corpus check` answers "can the next run actually succeed on this host?";
`corpus run --dry-run` prints the *plan* without writing artifacts. On a login node
with `$GROBID_URL` exported from step 6, these are the only lines that should not
be ✓ — two warnings and one informational:

| Line | Expected on the login node | Why |
|---|---|---|
| `GPU: none` | ⚠ | GPU phases run on `gpu_h200` / `gpu` via `sbatch`; the login node has no device |
| `Figure panels: would downgrade to OCR floor` | ⚠ | advisory only — `batch_pass3b.sh` re-decides on the GPU node, where `vision-local` is live |
| `OCR compression helpers: pngquant, jbig2enc absent` | • | informational; `ocrmypdf` drops `--optimize 2 → 1`, so output is larger, not wrong |
| everything else | ✓ | if not, stop and fix it before submitting |

In particular `Grobid`, `Taxonomy`, `bib`, `lexicon`, `Models` and the OCR toolchain
should all be ✓ by this point. Exit codes (#61): **0** green, **2** config error,
**3** precondition failure — so `corpus check` inside a `set -e` wrapper aborts the
submit. Since `check` honors `$GROBID_URL`, a ✗ Grobid here means Grobid really is
unreachable from this host, not that you're checking from the wrong node.

`corpus run --dry-run` ends with `reconcile failed in dry-run (exit 1) … Continuing`
on a corpuscle that hasn't been built yet. That is expected, not a preflight
failure — `reconcile` needs `biblio_authority.sqlite`, which only a real run
produces. The dry run still exits 0.

*(Steps 1–7 verified end-to-end on the cluster 2026-08-01, against
`corpuscles/siphonophore_20260731`.)*

## Smoke test: the gold corpuscle, before the full corpus

Before committing to the full-corpus run, smoke-test end-to-end on a small
sample. Unlike step 7's `corpus run --dry-run`, this is a real build: it submits
real CPU and GPU jobs and writes real artifacts, which is the point — you
inspect them below. It is also a different exercise from
[PLATFORM_SMOKE.md](PLATFORM_SMOKE.md), which checks that the *install* works
across platforms against the 4-paper demo corpus; this checks that *your*
production `config.yaml` and *your* PDFs survive the whole chain at small scale.

**The sample is the 35 transcribed documents**, not an ad-hoc slice. The
siphonophore library's `transcriptions/` tree holds verbatim, page-by-page
transcriptions of 35 documents (761 pages, 1594–2026, 13 languages), made from
rendered page images only and never by correcting an extractor's output. They
were chosen to span the collection's hardest axes — Fraktur, Cyrillic, CJK,
vertical Japanese, plate-only atlases, documents with no text layer at all —
which is exactly what a smoke test wants to exercise, and they come with ground
truth attached, so the same build serves three purposes: this rehearsal, the
release drift reference (#187), and extraction-accuracy scoring (#192).

**It is a cheaper build than the 30 PDFs it replaced, despite being 35** — 761
pages against 1,290, and ~644 pages needing OCR against ~916. The old sample was
never picked for speed: it carried a 235-page thesis plus four more 100-page
monographs, three of them full-page scans. `BATCH_SIZE` defaults to 64, so 35
PDFs is a single array task and wall-clock tracks total work.

The one heavy document is **Totton1965a** — 314 pages at 100% raster coverage,
so scan detection calls it a scan and `redo_ocr` re-OCRs every page whatever its
text layer claims. That is half the set's entire OCR load, and it is kept on
purpose: 314 scanned pages against 314 independently transcribed pages is the
largest OCR ground truth available, in the exact document class this corpuscle
exists to measure. **If it turns out to dominate wall-clock, split the roles
rather than dropping it** — keep all 35 for scoring and #187, run the rehearsal
on a subset. `input_pdfs` is a directory, so that costs a directory of symlinks
and no code change. Dropping it from *scoring* is the one thing not to do.

Under the config-driven flow (#138) it is just a second corpuscle — its own
directory + `config.yaml` pointing at a slice of PDFs — selected by exporting
`CORPUS_CONFIG`:

```bash
# 1. Materialise the 35 gold PDFs. sources.json maps each transcription stem
#    to its PDF and that file's sha256 — resolve through the file and verify
#    the checksum rather than guessing a shelf letter. (The library holds two
#    editions of Lery, both once named Lery1594.pdf; only the 1594 is
#    transcribed, and they are not interchangeable.)
mkdir -p "$BOUCHET_PROJECT/siphonophores_gold/library"
python - <<'PY'
import hashlib, json, pathlib, shutil
lib  = pathlib.Path("siphonophores")           # the library repo
dest = pathlib.Path("siphonophores_gold/library"); dest.mkdir(parents=True, exist_ok=True)
for stem, rec in json.loads((lib / "transcriptions/sources.json").read_text()).items():
    src = lib / rec["pdf"]
    got = hashlib.sha256(src.read_bytes()).hexdigest()
    assert got == rec["sha256"], f"{stem}: checksum mismatch, refusing to copy"
    shutil.copy2(src, dest / src.name)
PY

# 2. A gold corpuscle. Use today's date, e.g. 20260825.
mkdir -p "$BOUCHET_PROJECT/corpuscles/siphonophore_gold_YYYYMMDD"
cd "$BOUCHET_PROJECT/corpuscles/siphonophore_gold_YYYYMMDD"
# Copy the production config and repoint input_pdfs; output_dir: . keeps
# the sample's artifacts out of the production tree. `bib:` stays pointed at
# the production siphonophores.bib — per-document directives live there.
cp "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/config.yaml" .
#   edit:  input_pdfs: ../../siphonophores_gold/library
export CORPUS_CONFIG="$BOUCHET_PROJECT/corpuscles/siphonophore_gold_YYYYMMDD/config.yaml"

# 3. Run the phases, from the repo (the slurm/ paths are relative to it).
#    The simplest path is the orchestrator, which inherits $CORPUS_CONFIG
#    via --export=ALL and runs the whole chain:
cd "$BOUCHET_PROJECT/corpus"
bash slurm/batch_pipeline.sh

#    …or submit phases by hand (each reads $CORPUS_CONFIG):
# sbatch slurm/batch_process_corpus.sh     # extract (CPU, day)
# sbatch slurm/batch_pass3b.sh             # vision  (GPU, gpu_h200)
# sbatch slurm/batch_embed.sh              # embed   (GPU)
# sbatch slurm/batch_finalize.sh           # post + bundle (CPU)
```

Inspect `$BOUCHET_PROJECT/corpuscles/siphonophore_gold_YYYYMMDD/documents/<HASH>/summary.json` and `figures_report.html` for a handful of papers. Confirm:

- Grobid metadata populated (`metadata.json`)
- Figures extracted + captioned, panels listed in `figures.json`
- Pass 3b ROIs populated, compounds detected where expected
- Pass 3c renamed compound images to range notation (e.g., `fig_3-4.png`)
- `chunks.json` has reasonable content, embedded into `vector_db/lancedb`

**Expect worse rates here than on the production corpus, and don't loosen the
shipped ceilings for it.** `tests/test_corpus_wide.py`'s soft-rate ceilings are
calibrated across two ordinary corpora; a hardest-cases set will legitimately
exceed some of them. Use the `CORPUS_SOFT_RATE_CEILINGS` env override, which
exists for precisely this third corpus.

Only then submit the production run against the full input.

---

## Production run

Everything above is setup and rehearsal. This section is the real build, and it
is self-contained: with steps 1–7 done once, this is the only part you repeat.

### Environment

Almost everything is derived by [slurm/bouchet_paths.sh](../slurm/bouchet_paths.sh),
which every batch script sources. What you actually have to do:

```bash
module load miniconda
conda activate corpus          # required — see below
cd "$BOUCHET_PROJECT/corpus"   # the slurm/ paths are relative to the repo
```

**`conda activate corpus` is not optional, and not for the reason it looks.**
The orchestrator itself needs no env. But every phase script sources
`bouchet_paths.sh` *before* it runs `conda activate`, and that file sets
`LD_LIBRARY_PATH` from `$CONDA_PREFIX` — which at that moment is whatever the
*submitting* shell had, propagated by `sbatch --export=ALL`. Submit without
activating and the jobs run with no `LD_LIBRARY_PATH`, which is
[#20](https://github.com/caseywdunn/corpus/issues/20): docling and torch fall
back to system libraries and can produce mis-structured output without
crashing. `batch_process_corpus.sh` has a fail-loud import preflight that
catches the worst of it, but don't rely on it.

| Variable | Who sets it | Notes |
|---|---|---|
| `BOUCHET_PROJECT` | you, **only if** not `/nfs/roberts/project/pi_cwd7/cwd7` | everything else derives from it |
| `CORPUS_CONFIG` | you, **only** for a one-off build that isn't `corpuscles/current` | otherwise resolved from the symlink; hard error if neither exists |
| `HF_HOME` | `bouchet_paths.sh` → `$BOUCHET_PROJECT/cache/huggingface` | **don't set it by hand for the run.** You *do* need it exported for the manual `corpus prefetch` in §5, and it must match this path or the jobs re-download (the "Stale HF_HOME" pitfall) |
| `LD_LIBRARY_PATH` | `bouchet_paths.sh`, from `$CONDA_PREFIX` | see the warning above |
| `REPO_DIR`, `CACHE_DIR` | `bouchet_paths.sh` | derived from `BOUCHET_PROJECT` |
| `GROBID_URL` | `batch_pipeline.sh` | discovered from its own Grobid job and injected into extract. Any value you export is overwritten |

Optional knobs:

| Variable | Default | Effect |
|---|---|---|
| `NUM_BATCHES` | *auto* — `ceil(PDF files / BATCH_SIZE)` | Stage 1 CPU array tasks. Derived from the corpuscle's own `input_pdfs`, so you no longer have to set it; export a value only to override |
| `BATCH_SIZE` | `64` | PDFs per Stage 1 task — see "Why `BATCH_SIZE` is 64" below |
| `NUM_PASS3B_BATCHES` | `1` | Pass 3b GPU array tasks. Deliberately *not* tied to `NUM_BATCHES` |
| `PASS3B_BATCH_SIZE` | `256` | papers per Pass 3b task |
| `HF_HUB_OFFLINE` | unset | `1` pins the run to the cached model snapshot; makes a build reproducible against a moving upstream (§5) |
| `ENRICH_BHL` | unset | `1` adds BHL enrichment in finalize — slow, rate-limited, many hours |
| `SKIP_BUNDLE` | unset | `1` runs the cross-paper DBs without distilling `_serve/` |

### Before you launch

**Run it from the login node.** `batch_pipeline.sh` submits seven jobs, polls
`squeue` and `curl` while sleeping, prints a summary and exits — it does no
compute, and every phase runs in its own SLURM job with `afterok` dependencies,
so nothing depends on the launcher staying alive. Under `salloc` you would hold
an allocation to run a poller, and if that allocation expires during the ~15 min
Grobid wait you get an orphaned Grobid with no cleanup job attached to it.

It is also **not resumable**: if it dies partway, some jobs are submitted and
some aren't, and re-running submits a second Grobid plus a duplicate chain.
If that happens, `squeue --me`, cancel everything it created, and start over.

**Cancel any Grobid you started by hand.** `slurm/batch_pipeline.sh` submits its
*own* Grobid unconditionally — there is no way to hand it an existing one. It
overwrites `$GROBID_URL` with its own job's node, and the cleanup job it
schedules cancels only the job it started. A Grobid left over from §6 or §7 will
therefore sit idle holding a `week` allocation for the full 48 hours.

```bash
# Find it — the job name is `grobid`:
squeue -u "$USER" -o "%.10i %.12j %.8T %.10M %N"

scancel <grobid_job_id>
```

If you kept the shell from §6, `scancel "$GROBID_JOB"` does it directly.

**Confirm the thing that still fails silently rather than loudly:**

```bash
# A leftover `export CORPUS_CONFIG=…gold…` from a smoke test will
# redirect a production submit into the sample corpuscle.
echo "${CORPUS_CONFIG:-<unset — will use corpuscles/current>}"
ls -l "$BOUCHET_PROJECT/corpuscles/current"
```

Under-provisioning `NUM_BATCHES` used to belong on this list — it left the
tail of the library unprocessed with no error. It no longer can: the
launcher derives the array size from the corpuscle's `input_pdfs` and
prints what it decided, e.g.

```
Auto-sized Stage 1: <N> PDFs / 64 per task = <ceil(N / 64)> tasks
```

Check that line against the library size if you want the belt-and-braces
version; the count comes from the same `find` as the header of this doc.

### Launching

The pipeline orchestrator (`slurm/batch_pipeline.sh`) handles Grobid startup, extract job-array submission, Grobid cleanup, and vision + embed + finalize chaining automatically — every job runs `corpus run --only <phase>` against `$CORPUS_CONFIG`.

The orchestrator parallelizes both Stage 1 (CPU) and Pass 3b (GPU)
independently. Pass 3b looks like the long pole and is not: only figures
whose caption declares more than one panel reach the VLM, which on one
reference build was about 900 of ~22,000 figure records (4.3%).

That reference run took 5 h 10 m end to end. These are illustrative
measurements; task count is `ceil(unique PDFs / 64)`, and elapsed time scales
with the current paper, chunk, and eligible-figure counts:

| Phase | Shape | Elapsed |
|---|---|---|
| extract | `ceil(unique PDFs / 64)` tasks | 2 h 35 m (slowest task in the reference run; most finished in 1–4 min) |
| vision (Pass 3b) | 1 GPU task | 1 h 24 m |
| embed | 1 GPU task | 21 m |
| finalize | 1 CPU task | 1 h 03 m |

Extract's wallclock is set by whichever slice draws the scanned Fraktur and
Cyrillic monographs, not by the average paper. Pass 3b's is dominated by the
walk over every document directory rather than by the 934 VLM calls, so it
scales with corpus size even though the GPU work does not.

The parallelism knobs are in the [Environment](#environment) table above.

```bash
cd "$BOUCHET_PROJECT/corpus"

# The whole build, hands-off. NUM_BATCHES is derived from the corpuscle's
# own input_pdfs, so there is nothing to recompute as the library grows.
bash slurm/batch_pipeline.sh

# Override the slice size (NUM_BATCHES re-derives from it automatically):
BATCH_SIZE=128 bash slurm/batch_pipeline.sh

# Fan Pass 3b out across GPUs. Rarely needed — only for a genuinely
# figure-bound corpus. Mind the 16-GPU per-user cap and check
# `sinfo -p gpu_h200` first.
NUM_PASS3B_BATCHES=4 bash slurm/batch_pipeline.sh
```

### Why `BATCH_SIZE` is 64 and not 256

The 2026-08-02 siphonophore build ran 8 tasks × 256 PDFs against Stage 1's
24 h wall. Two tasks drew the OCR-heavy scans and hit the wall at 252/256
and 225/256; the tasks that *did* finish cleared it by as little as 55 min.
Because Pass 3b, Embed and Finalize all chain on `afterok`, those two
TIMEOUTs cancelled the entire downstream chain — silently, with no log, no
mail and no queue entry. See "Chain integrity" below.

64 PDFs/task gives roughly a 3× margin at the worst observed rate, and more
tasks spread the OCR-heavy tail instead of concentrating it in one slice.
Empty and already-complete tasks cost ~4 min each, so over-provisioning is
cheap; under-provisioning silently leaves papers unprocessed.

### The OCR worker pool is sized from the allocation, not the node

Bouchet's stage1 nodes are large — `a1130u24n01` reports `CPUTot=64` and
991 GiB — and Stage 1 requests `--cpus-per-task=8`. ocrmypdf, left to
itself, reads `multiprocessing.cpu_count()` and sees the 64. It then runs
64 Tesseract workers inside an 8-CPU cgroup, each page gets an eighth of a
core, and pages cross `ocr.tesseract_page_timeout` — at which point
ocrmypdf copies the un-OCR'd image into the output and exits 0. The page
survives visually with an empty text layer.

Two full builds of the same library lost text this way on ~9.5% of
documents each. It is load-dependent, so it takes a *different* set every
run: 31 documents lost >80% of their text between one build and the next
while 28 independently gained it back. Nothing OOMs — 64 × 2.5 GB fits
comfortably inside a 256 G request — which is why a memory-based cap never
fired (#254).

`corpus` now reads `$SLURM_CPUS_PER_TASK` (then the affinity mask, then a
cgroup quota) and passes `--jobs` explicitly, so nothing is needed in the
config for this. Check the `Running OCR on …` line in a Stage 1 log: it
prints `jobs=` and it should equal `--cpus-per-task`. Pin `ocr.jobs` in
the corpuscle config only if you are running an older build, or if you
change `--cpus-per-task` in a way the env doesn't reflect.

Blanked pages are now recorded per document as
`scan_detection.json`'s `pages_blanked` and gated at `error` as
`ocr_pages_blanked`, so a build that still hits this says so:

```bash
corpus status --filter-gate ocr_pages_blanked --list-hashes
corpus run --only extract --re-process-flagged ocr_pages_blanked
```

### Chain integrity

`afterok` on a job *array* requires **every** task to exit 0. One TIMEOUT
takes out Pass 3b, Embed and Finalize. Three things now guard this:

- `NUM_BATCHES` is auto-derived, so the array can no longer be too small
  for the library.
- Stage 1 carries `--mail-type=FAIL,TIME_LIMIT`, so a wall hit sends mail.
- A `chain-watchdog` job runs on `afterany` after Stage 1 and always
  reports: the per-task states, the completed-document count, and an
  explicit warning naming the cancelled downstream job IDs. Read
  `logs/slurm-watchdog-<jobid>.out` before assuming a build finished.

`batch_pipeline.sh` also prints `squeue --me -o "%.12i %.18j %.10T %.30E"`
right after submit — every downstream job should show a pending
`Dependency` reason there.

The orchestrator:
1. Starts Grobid as a SLURM job, waits for it to be alive, exports `$GROBID_URL` for the extract job
2. Submits **extract** (`batch_process_corpus.sh`) as a job array (`--array=0-N`), each task processing a deterministic slice of the sorted hash list
3. Schedules Grobid cleanup after extract completes (runs regardless of success/failure)
4. Queues **vision** (`batch_pass3b.sh`, GPU) and **embed** (`batch_embed.sh`, GPU) with `afterok` dependency on the full extract array
5. Queues **finalize** (`batch_finalize.sh` — `corpus run --only post` then `--only bundle`) with `afterok` on **both** embed and vision. Both are needed: the served bundle copies `figures.json` and `figures/*.png`, which Pass 3b rewrites and Pass 3c renames, so gating on embed alone let `bundle` capture pre-vision ROIs and stale figure filenames
6. Queues a **chain-watchdog** on `afterany` after extract, which reports the array outcome whether or not the chain survived

So the orchestrator now runs the whole build, including the cross-paper DBs and served-bundle distill, hands-off. Every job invokes `corpus -c "$CORPUS_CONFIG" run --only <phase>`.

For manual submission without the orchestrator (each phase reads `$CORPUS_CONFIG`; resume is implicit):

```bash
cd "$BOUCHET_PROJECT/corpus"
export GROBID_URL=http://<grobid_node>:8070       # extract needs Grobid (step 6)
S1=$(sbatch --parsable --array=0-27 slurm/batch_process_corpus.sh)   # extract
P=$(sbatch --parsable --dependency=afterok:$S1 slurm/batch_pass3b.sh)  # vision
E=$(sbatch --parsable --dependency=afterok:$S1 slurm/batch_embed.sh)   # embed
sbatch --dependency=afterok:$E:$P slurm/batch_finalize.sh              # post + bundle → _serve/
```

Size `--array` as `ceil(PDF files / BATCH_SIZE) - 1`; the orchestrator does this for you.

Resume is implicit in `corpus run`, so restarts are cheap — re-queuing a phase re-processes only papers whose inputs changed. Without `NUM_PASS3B_BATCHES` (or set to 1), Pass 3b runs as a single job, which is right for all but genuinely figure-bound corpora.

**Do not trust `corpus run --only embed --dry-run`.** The real run gates resume on the marker's `embedding_model` + `embedding_dim` (`_marker_matches_backend`, `pipeline/embed.py:216`), but the dry-run path checks only that `vector_db/<HASH>_embedded.done` exists (`embed.py:296-306`). Stage 1 writes a *chunking* receipt at that same path (the `ingest_to_vector_db` stub, `pipeline/chunking.py:149`), which carries neither key. So on a corpus that has finished extract but never embedded, the dry-run reports "already have a marker" for every paper while the real run correctly embeds them all. Check `vector_db/lancedb/` for the actual state — if that directory is absent, nothing has been embedded regardless of how many `.done` files there are.

## Post-pipeline cross-paper databases + served bundle

The extract/vision/embed phases produce per-paper artifacts under
`<output_dir>/documents/<HASH>/`. The **post** phase then runs the four
cross-paper builds in dependency order, and the **bundle** phase distills
the served bundle. Both are wrapped by `slurm/batch_finalize.sh`, which the
orchestrator chains automatically after embed — so normally you do nothing
here. The four builds (no longer separate top-level scripts; #138 folds
them into `corpus run --only post`) are:

1. **biblio authority** — deduplicated works + citation graph (`bib.authority`)
2. **taxon mentions** — cross-paper taxon-mention index (`pipeline.taxon_mentions`)
3. **in-text citations** — TEI body → `intext_citations.json` per paper (`pipeline.intext_citations`)
4. **reconcile** — merge ghost cited-references onto corpus papers (`bib.reconcile`)

To run the tail by hand (e.g. after fixing a build issue):

```bash
cd "$BOUCHET_PROJECT/corpus"
export CORPUS_CONFIG="$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/config.yaml"

# Cross-paper DBs only. ENRICH_BHL=1 adds Biodiversity Heritage Library
# coverage for pre-DOI refs (slow, rate-limited ~1 req/s — many hours on
# the week partition); omit for a fast build without BHL.
sbatch slurm/batch_finalize.sh
# ENRICH_BHL=1 sbatch slurm/batch_finalize.sh

# Cross-paper DBs but skip the bundle distill:
# SKIP_BUNDLE=1 sbatch slurm/batch_finalize.sh
```

`batch_finalize.sh` runs `corpus run --only post` then `corpus run --only
bundle`. The post phase honors `--force-rebuild*` only when you ask; a plain
re-run is idempotent. Confirm the cross-paper SQLites landed:

```bash
ls -la "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD"/{taxonomy,biblio_authority,taxon_mentions}.sqlite
```

The bundle phase distills into `<output_dir>/_serve/` — for the production
corpuscle that's `$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/_serve/`. This
replaces the old standalone `package_for_serve.py` / `SERVE_BUNDLE_DIR`
step: `_serve/` **is** the deployable artifact (path-scrubbed + audited +
manifested). Re-distill any time with:

```bash
corpus -c "$CORPUS_CONFIG" run --only bundle
cat "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD/_serve/bundle_manifest.json"
```

`_serve/` is what gets uploaded to S3 and consumed by the EC2 deploy — see [DEPLOY.md](../DEPLOY.md) §5 (the host pulls it with `deploy/update.sh <version>`); §6 is the post-deploy smoke test.

## Acceptance testing a completed build

After `batch_finalize.sh` completes, verify the corpuscle on an **interactive compute
node** (not the login node — the MCP server loads the ~600 MB BGE-M3 embedder when
the first semantic-search tool is called, which requires more memory than the shared
login node allows):

```bash
salloc -p devel -t 2:00:00 --mem=32G
module load miniconda
conda activate corpus

# 1. Programmatic smoke test (corpus-agnostic, any corpuscle)
python "$BOUCHET_PROJECT/corpus/tools/smoke_test_sse.py" \
    "$BOUCHET_PROJECT/corpuscles/siphonophore_YYYYMMDD" \
    --server-startup-timeout 120

# 2. Manual acceptance prompts (siphonophore-specific — see note)
#    Connect Claude to the running MCP server and paste each prompt from:
#    dev_docs/ACCEPTANCE_PROMPTS.md
```

> **Note:** `ACCEPTANCE_PROMPTS.md` is a siphonophore-specific example.
> For a different corpus, copy it and substitute taxon/author names.

## Partition reference

| Phase (`--only`) | Script | Partition | GPU? | Walltime |
|---|---|---|---|---|
| `extract` (OCR + docling + Grobid + Pass 2.5) | `slurm/batch_process_corpus.sh` | `day` | no | 24 h |
| `vision` (Pass 3b + 3c, Qwen2.5-VL-7B) | `slurm/batch_pass3b.sh` | `gpu_h200` | 1 | 24 h |
| `embed` (BGE-M3) | `slurm/batch_embed.sh` | `gpu` | 1 | 4 h |
| `post` + `bundle` (cross-paper DBs + served bundle) | `slurm/batch_finalize.sh` | `day` | no | 12 h |
| (Grobid service) | `slurm/batch_grobid.sh` | `week` | no | 48 h |

Each script runs `corpus -c "$CORPUS_CONFIG" run --only <phase>`. Adjust walltimes as the corpus grows.

Grobid sits on `week`/48 h rather than `day`/24 h **on purpose**: it is submitted
before extract, so an equal wall guaranteed it died first. Papers extract handles
after Grobid dies get placeholder metadata, and implicit resume will not retry them
(their inputs are unchanged) — a silent quality loss. It costs nothing to outlive
extract, because the `afterany` grobid-cancel job tears the service down the moment
extract ends. **Keep Grobid's walltime strictly greater than extract's.**

The per-user QoS caps that actually bind these submissions (`sacctmgr show qos`,
2026-08-01):

| Partition | Max wall | Per-user cap |
|---|---|---|
| `day` | 1 d | 1000 CPUs, 15 TiB mem |
| `week` | 7 d | 96 CPUs, 1.5 TiB mem |
| `devel` | 6 h | 4 CPUs, 60 GiB |
| `gpu`, `gpu_h200` | 2 d | 16 GPUs |
| `gpu_devel` | 6 h | 2 GPUs, 12 CPUs, 256 GiB |

Both GPU phases sit well inside the 16-GPU cap at the defaults, which run
Pass 3b as a single job. The 16-GPU cap is why `NUM_PASS3B_BATCHES` no
longer inherits `NUM_BATCHES`: at `BATCH_SIZE=64` the extract array is
~28 tasks, and matching that on `gpu_h200` would queue against the cap for
work one GPU finishes in about an hour. The extract array itself is bounded
by `day`'s 1000-CPU cap (28 × 8 = 224 CPUs), not by anything in the scripts.

Note that the `gpu` partition is **heterogeneous** — it carries `a40` (48 GB),
`l40s` (48 GB) and `rtx_5000_ada` (32 GB), so `batch_embed.sh`'s bare
`--gpus=1` can land on any of the three. BGE-M3 fits comfortably on all of them.
`gpu_devel` carries one node each of `h200`, `b200`,
`rtx_pro_6000_blackwell` and `rtx_5000_ada`, which is why it is the right
partition for interactive GPU checks.

## Common pitfalls

- **`jobs=` in the Stage 1 log should equal `--cpus-per-task`.** If it says 64 on an 8-CPU task, the allocation isn't being read and OCR is oversubscribing the cgroup — pages will be silently blanked rather than fail. See "The OCR worker pool is sized from the allocation" above.
- **Missing Tesseract language packs.** The most likely way to get a subtly bad build. `conda env create` alone leaves you with **English-only** OCR; `bash tools/install_tessdata.sh` is a required setup step (see §2). Non-English papers don't error — they just OCR badly as English. Check with `ls $CONDA_PREFIX/share/tessdata/`.
- **Stale `HF_HOME`.** If a job re-downloads a model, `HF_HOME` isn't being honored — check that the export in the SLURM script points to a path you actually populated.
- **Grobid URL.** SLURM compute nodes can't talk to your laptop's `localhost:8070` — `$GROBID_URL` (which overrides the config's `grobid.url`) must resolve to a host visible from the job's node. If the Grobid node goes down mid-run, subsequent papers get placeholder metadata; a re-run's implicit resume won't retry them unless their inputs changed — force it with `corpus run --only extract --re-process-flagged <gate>` or by deleting the affected `metadata.json`.
- **`corpus check` and Grobid.** `check` honors `$GROBID_URL` exactly as `run` does, so the two always agree about which Grobid they mean. If `check` reports `localhost:8070`, `GROBID_URL` simply isn't exported in that shell — the check is right, and the extract job you were about to submit would have used the same wrong address. Re-do the export from §6 and re-run.
- **Config not found / wrong corpuscle.** Every phase script reads `$CORPUS_CONFIG`. If a job dies with "no config.yaml" or builds the wrong tree, confirm `$CORPUS_CONFIG` points at the intended `config.yaml` (default is the production corpuscle; a leftover `export CORPUS_CONFIG=…gold…` from a smoke test will silently redirect a production submit).
- **LFS on extract.** If extract can't see the full PDFs (only LFS pointers), re-run `git lfs pull` in `$BOUCHET_PROJECT/siphonophores`.
- **Don't downgrade torch to "fix" a GPU problem.** The pinned `torch==2.12.0` works on every GPU type on the cluster — all of them run driver 580.159.04 / CUDA 13.0. Older cu128 builds cannot target the Blackwell cards (`b200` sm_100, `rtx_pro_6000_blackwell` sm_120) at all, so a well-meant downgrade silently costs you that hardware and breaks the #98 pin. `batch_pass3b.sh` targets `gpu_h200` for VRAM headroom and throughput on Qwen2.5-VL-7B, **not** because other cards fail; its preflight `torch.cuda.is_available()` guard stays as a cheap regression check.
- **lmod module cache flakiness.** `module load CUDA/12.6.0` occasionally errors with `CUDA/12.6.0.lua: Empty or non-existent file` even though `module spider CUDA` lists it. `slurm/batch_pass3b.sh` skips the explicit CUDA module entirely — torch ships bundled CUDA userspace libs in `site-packages/nvidia/`, so no CUDA module is needed for the GPU phases at all. If another script hits this, retry with `module --ignore-cache load …`.

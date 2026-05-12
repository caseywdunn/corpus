# Platform-portability smoke test

Pre-release sanity check that the supported matrix in
[README.md](../README.md#supported-platforms) actually works
end-to-end. Run before tagging a new release.

The test is the same on every target: build the demo corpus with the
pipeline + bundle distillation, then verify `bundle_info` returns
expected values from the served bundle. Differences between targets
live in the install path (conda env, system tools) and the Grobid
network topology, not in the success criteria.

## Targets

| # | Target | Where |
|---|---|---|
| 1 | macOS arm64 | MacBook bare metal |
| 2 | linux-x86_64 | Bouchet (or any Linux box with Docker) |
| 3 | linux-arm64 | Docker on the MacBook, native (no emulation) |
| 4 | linux-x86_64 (clean-room) | Docker on the MacBook under Rosetta — optional, validates env file from a clean container |

macOS x86_64 is explicitly unsupported — no torch ≥ 2.4 wheels exist
for that arch. Don't smoke-test it.

linux-aarch64 outside Docker (e.g. a Graviton EC2 host) is not in this
runbook today; #3 covers the architecture but in a different distro.
Promote later if a real Graviton target appears.

## Success criteria

Each target must:

1. `conda env create -f environment.yaml` (or `docker build`) finishes
   without error.
2. `pip install -e .` finishes without error.
3. `corpus check` reports Grobid reachable + config valid.
4. `corpus run` on the bundled `demo/` corpus completes — all 11 PDFs
   reach `pipeline_state.json` status `done` (`corpus status --report`
   shows 11 / 11).
5. `bundle_manifest.json` is written under `demo/output/_serve/`,
   contains `paper_count: 11`, and the absolute-path audit logs
   `Path scrub: rewrote N files; audit clean.` (covers
   [#70](https://github.com/caseywdunn/corpus/issues/70)).
6. `corpus serve --output-dir demo/output/_serve` starts, and
   `bundle_info` via any MCP client returns the same `paper_count` +
   the bundle version stamped in `pipeline/version.py`.

If any of those fails, the release is blocked.

The vision pass is `--no-vision` everywhere in this runbook — vision
backends are GPU-bound and not testable from Mac Docker (no GPU
passthrough on Docker Desktop). Vision-pass quality is validated
separately on Bouchet.

## (1) macOS arm64 — bare metal

Prerequisites already met if INSTALL.md was followed once. Sanity-
check arch and run:

```bash
~/miniforge3/envs/corpus/bin/python -c "import platform; print(platform.machine())"
# expect: arm64

docker compose up -d grobid                # linux/amd64 image, runs under Rosetta
curl -fsS http://localhost:8070/api/isalive  # expect: true

cd demo && corpus run --no-vision
corpus status --report
corpus serve --output-dir output/_serve &
# point an MCP client at it, call bundle_info, kill the server
```

## (2) linux-x86_64 — Bouchet

Same conda env, same demo, on the production HPC target. The SLURM
chain isn't part of the smoke test (separately covered by BOUCHET.md);
this is the laptop-equivalent path on a Linux box.

```bash
module load miniconda                      # YCRC convention
conda env create -f environment.yaml       # or `conda env update -f environment.yaml --prune` to refresh
conda activate corpus
pip install -e .
bash tools/install_tessdata.sh

# Grobid via Apptainer/Singularity — see INSTALL.md#grobid-on-bouchet
singularity build grobid.sif docker://lfoppiano/grobid:0.8.1
singularity run --bind $HOME grobid.sif &
curl -fsS http://localhost:8070/api/isalive

cd demo && corpus run --no-vision
corpus status --report
```

On any non-HPC Linux box, swap the apptainer line for
`docker compose up -d grobid`.

## (3) linux-arm64 — Docker on Apple Silicon (native)

The fastest of the three Docker variants because the container's arch
matches the host's. Use this as the routine pre-release check for
linux-arm64.

```bash
# One-time: build the smoke image natively for linux/arm64.
docker buildx build --platform=linux/arm64 \
    -f dev_docs/Dockerfile.smoke \
    -t corpus-smoke:arm64 --load .

# Grobid runs as a separate linux/amd64 container under Rosetta; the
# corpus container reaches it via host.docker.internal because
# Docker Desktop's bridge network isolates container localhost.
docker compose up -d grobid
curl -fsS http://localhost:8070/api/isalive

# Run the smoke. --grobid-url overrides demo/config.yaml's
# http://localhost:8070 (which would be the container's own loopback)
# with the host gateway alias.
docker run --rm -it \
    -v "$PWD":/corpus \
    -w /corpus/demo \
    corpus-smoke:arm64 \
    bash -c "corpus run --no-vision --grobid-url http://host.docker.internal:8070 \
             && corpus status --report"

# Inspect bundle outside the container — demo/output/ lives on the host
# bind mount.
jq '.paper_count' demo/output/_serve/bundle_manifest.json   # expect: 11
```

## (4) linux-x86_64 — Docker on Apple Silicon (Rosetta, optional)

Same Dockerfile, x86_64 platform. Slower than (3) (emulation), but
validates that the environment.yaml install path works on
linux-x86_64 from a clean container — independent of any state on the
long-lived Bouchet env.

```bash
docker buildx build --platform=linux/amd64 \
    -f dev_docs/Dockerfile.smoke \
    -t corpus-smoke:amd64 --load .

# Grobid + the smoke command are identical to (3); only the image tag
# differs.
docker run --rm -it \
    -v "$PWD":/corpus \
    -w /corpus/demo \
    corpus-smoke:amd64 \
    bash -c "corpus run --no-vision --grobid-url http://host.docker.internal:8070 \
             && corpus status --report"
```

## Troubleshooting

**`Dynamo is not supported on Python 3.12+`** during model load on
macOS — the active conda is x86_64 (Intel anaconda under Rosetta), not
arm64. Recreate the env with miniforge. See
[INSTALL.md](../INSTALL.md#apple-silicon-use-miniforge-not-intel-anaconda).

**`pngquant not on PATH`** warning during OCR — expected on macOS unless
you `brew install pngquant`. Pipeline auto-degrades `--optimize` from
2 to 1; not a smoke-test failure on its own.

**`OSError: [Errno 14] Bad address: 'g++'`** in a docling crash — only
surfaces if `TORCH_COMPILE_DISABLE=1` got unset somewhere. The default
in `pipeline/__init__.py` covers this; check whether the host shell
has explicitly set `TORCH_COMPILE_DISABLE=0`.

**`Table 'document_chunks' already exists`** on a re-run against an
existing `output/` — the existing `pipeline/embed.py:312` guards this
with an open-if-exists pattern. If you hit it, you're on an old clone
(pre-v0.3.0); pull and retry. Tracked at
[#71](https://github.com/caseywdunn/corpus/issues/71).

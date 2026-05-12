# Platform-portability smoke test

Pre-release sanity check that the supported matrix in
[README.md](../README.md#supported-platforms) actually works
end-to-end. Run before tagging a new release. Both targets must pass;
either failing blocks the release.

The test is the same on every target: build the demo corpus with the
pipeline + bundle distillation, then verify `bundle_info` returns
expected values from the served bundle. Differences between targets
live in the install path (conda env, system tools, Grobid topology),
not in the success criteria.

## Targets

| # | Target | Where |
|---|---|---|
| 1 | macOS arm64 | MacBook bare metal (miniforge corpus env) |
| 2 | linux-x86_64 | Bouchet — **clean env recreate**, not the long-lived one |

macOS x86_64 is explicitly unsupported — no torch ≥ 2.4 wheels exist
for that arch. linux-aarch64 is not currently supported either; add a
target here when it does.

Step (2) intentionally tears down the existing Bouchet conda env so
the test catches `environment.yaml` regressions instead of riding on
months of accumulated state. If you'd rather keep your working env
intact, swap `-n corpus` for `-n corpus-smoke` in the recreate step —
you still get the clean-from-environment.yaml signal without nuking
the live env.

**Highest-signal alternative for target (2):** run
[`dev_docs/ec2_smoke.sh`](ec2_smoke.sh) on a clean Ubuntu EC2
instance. Same install path, but from absolutely nothing — no
pre-existing conda, Docker, or HF cache to ride on. The script also
emits a programmatic pass/fail summary against every criterion below,
so the release operator gets a single `exit 0` vs. `exit 1` instead
of having to eyeball logs. ~20–25 min wall time end-to-end, ~$2–3
EC2 cost.

## Success criteria

Each target must:

1. `conda env create -f environment.yaml` finishes without error from
   a clean state (no pre-existing `corpus` env).
2. `pip install -e .` finishes without error.
3. `corpus check` reports Grobid reachable + config valid.
4. `corpus run --no-vision` on the bundled `demo/` corpus completes —
   all 11 PDFs reach `pipeline_state.json` status `done`
   (`corpus status --report` shows 11 / 11).
5. `bundle_manifest.json` is written under `demo/output/_serve/`,
   contains `paper_count: 11`, and the absolute-path audit logs
   `Path scrub: rewrote N files; audit clean.` (covers
   [#70](https://github.com/caseywdunn/corpus/issues/70)).
6. `corpus serve --output-dir demo/output/_serve` starts, and
   `bundle_info` via any MCP client returns the same `paper_count` +
   the bundle version stamped in `pipeline/version.py`.

The vision pass is `--no-vision` everywhere here — vision backends
are GPU-bound and validated separately on Bouchet's `gpu_h200`
partition. The CPU-portability check is what this runbook owns.

## (1) macOS arm64 — bare metal

```bash
# Verify the active env is actually arm64. If this prints x86_64,
# stop and rebuild with miniforge — see
# INSTALL.md#apple-silicon-use-miniforge-not-intel-anaconda.
~/miniforge3/envs/corpus/bin/python -c "import platform; print(platform.machine())"
# expect: arm64

# Clean env recreate. Drop --force if you'd rather be prompted.
conda env remove -n corpus --yes
conda env create -f environment.yaml
conda activate corpus
pip install -e .
bash tools/install_tessdata.sh

docker compose up -d grobid                # linux/amd64 image, Rosetta
curl -fsS http://localhost:8070/api/isalive  # expect: true

cd demo && corpus -v check                 # -v required to see the ok lines
corpus -v run --no-vision                  # ~25–30 min total wall time on
                                           # an M-series MacBook; the WoRMS
                                           # taxonomy ingest is the long
                                           # pole (~10 min), then extract
                                           # (~6 min) + embed (~30s) + bundle.
corpus status --report                     # expect: 11 / 11 done
jq '.paper_count' output/_serve/bundle_manifest.json   # expect: 11

# Round-trip the MCP bundle_info tool against a freshly-served bundle.
# tools/smoke_test_sse.py spawns its own server on the requested port,
# initializes the MCP client, and calls bundle_info + list_papers.
python tools/smoke_test_sse.py demo/output/_serve --port 18080
# expect: "All layers passed." with bundle_version matching pipeline/version.py
```

## (2) linux-x86_64 — Bouchet (clean env)

```bash
module load miniconda                      # YCRC convention

# Clean env recreate — the whole point of this leg.
conda env remove -n corpus --yes
conda env create -f environment.yaml
conda activate corpus
pip install -e .
bash tools/install_tessdata.sh

# Grobid via Apptainer/Singularity — see INSTALL.md#grobid-on-bouchet
singularity build --force grobid.sif docker://lfoppiano/grobid:0.8.1
singularity run --bind $HOME grobid.sif &
curl -fsS http://localhost:8070/api/isalive

cd demo && corpus -v check                 # -v required to see the ok lines
corpus -v run --no-vision                  # wall time depends on Bouchet
                                           # load + WoRMS API rate; budget
                                           # 30–45 min for the demo.
corpus status --report                     # expect: 11 / 11 done
jq '.paper_count' output/_serve/bundle_manifest.json   # expect: 11
python tools/smoke_test_sse.py demo/output/_serve --port 18080
# expect: "All layers passed."
```

The SLURM chain (`slurm/batch_pipeline.sh`) isn't part of this
runbook — it's covered separately by
[BOUCHET.md](BOUCHET.md). The interactive `corpus run` here exercises
the same code path on a single login-node session and is fast enough
for a pre-release gate.

## Troubleshooting

**`Dynamo is not supported on Python 3.12+`** during model load on
macOS — the active conda is x86_64 (Intel anaconda under Rosetta),
not arm64. Recreate the env with miniforge. See
[INSTALL.md](../INSTALL.md#apple-silicon-use-miniforge-not-intel-anaconda).

**`pngquant not on PATH`** warning during OCR — expected on macOS
unless you `brew install pngquant`. Pipeline auto-degrades
`--optimize` from 2 to 1; not a smoke-test failure on its own.

**`OSError: [Errno 14] Bad address: 'g++'`** in a docling crash —
only surfaces if `TORCH_COMPILE_DISABLE=1` got unset somewhere. The
default in `pipeline/__init__.py` covers this; check whether the host
shell has explicitly set `TORCH_COMPILE_DISABLE=0`.

**`Table 'document_chunks' already exists`** on a re-run against an
existing `output/` — guarded by an open-if-exists pattern in
`pipeline/embed.py`. If you hit it, you're on an old clone (pre-v0.3.0);
pull and retry. Tracked at
[#71](https://github.com/caseywdunn/corpus/issues/71).

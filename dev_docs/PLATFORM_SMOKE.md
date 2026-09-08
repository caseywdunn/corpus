# Platform-portability smoke test — manual fallback / release-time verification

> **Authoritative coverage now lives in GitHub Actions.** As of #75, T1
> (Linux + Grobid) and T2 (macOS arm64) in
> [`.github/workflows/integration.yml`](../.github/workflows/integration.yml)
> exercise the same demo build on every push and every PR — that's
> where bundle-audit, manifest-shape, and SSE-round-trip regressions
> get caught now.
>
> This runbook is the manual fallback used at release time (when the CI
> tiers have already passed but the release operator wants a clean-env
> recreate signal that `actions/cache@v4` deliberately hides) and the
> escape hatch for paths that the GHA tiers can't reach: the EC2 clean-
> room (T3-bare — [`ec2_smoke.sh`](ec2_smoke.sh)) and bare-metal Bouchet.

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
4. `corpus run --no-vision` on the bundled 4-paper `demo/` corpus completes —
   all four PDFs reach `pipeline_state.json` status `done`
   (`corpus status --report` reports 4 documents and every stage row at
   `N / N`). How many stage rows there are is a property of the
   configuration, not of a healthy build: `--no-vision` records no
   `figure_pass*` rows, so do not assert a fixed row count.
5. `bundle_manifest.json` is written under `demo/output/corpus_bundle/`,
   contains `paper_count: 4`, and the absolute-path audit logs
   `Path scrub: rewrote N files; audit clean.` (covers
   [#70](https://github.com/caseywdunn/corpus/issues/70)).
6. `corpus serve --output-dir demo/output/corpus_bundle` starts, and
   `bundle_info` via any MCP client returns the same `paper_count` +
   the bundle version stamped in `pipeline/version.py`.

The vision pass is `--no-vision` everywhere here — vision backends
are GPU-bound and validated separately on Bouchet's `gpu_h200`
partition. The CPU-portability check is what this runbook owns.

## (1) macOS arm64 — bare metal

```bash
# Verify the active env is actually arm64. If this prints x86_64,
# stop and rebuild with miniforge — see
# INSTALL.md#apple-silicon-arm64-native-conda-required.
~/miniforge3/envs/corpus/bin/python -c "import platform; print(platform.machine())"
# expect: arm64

# Clean env recreate. Drop --force if you'd rather be prompted.
conda env remove -n corpus --yes
conda env create -f environment.yaml
conda activate corpus
pip install -e .
bash tools/install_tessdata.sh
bash tools/install_ocr_extras.sh

docker compose up -d grobid                # linux/amd64 image, Rosetta
# Startup is slower under Rosetta than native, so wait rather than probe once.
until [ "$(curl -fs http://localhost:8070/api/isalive 2>/dev/null)" = true ]; do
  sleep 5; echo "waiting for grobid..."
done                                       # expect: true, within ~1-2 min

cd demo && corpus -v check                 # -v required to see the ok lines
corpus -v run --no-vision                  # ~25–30 min total wall time on
                                           # an M-series MacBook; the WoRMS
                                           # taxonomy ingest is the long
                                           # pole (~10 min), then extract
                                           # (~6 min) + embed (~30s) + bundle.
corpus status --report                     # expect: 4 / 4 done
jq '.paper_count' output/corpus_bundle/bundle_manifest.json   # expect: 4

# Round-trip the MCP bundle_info tool against a freshly-served bundle.
# tools/smoke_test_sse.py spawns its own server on the requested port,
# initializes the MCP client, and calls bundle_info + list_papers.
python tools/smoke_test_sse.py demo/output/corpus_bundle --port 18080
# expect: "All layers passed." with bundle_version matching pipeline/version.py
```

## (1a) macOS arm64 — the local VLM's dtype (#258)

Open question, and the only thing in the tracker blocked on hardware rather
than on a decision. `figures.vision_dtype` defaults to `auto`, which is
float32 on MPS — **~30.4 GB of weights for Qwen2.5-VL-7B against ~15.2 GB at
half precision**, so it does not fit a 32 GB machine at all today. Apple
Silicon supports bf16 and fp16; what is unknown is whether Qwen2.5-VL is
numerically sound in either on Metal, and which of the two this model prefers.
Nobody has run it.

The default does not move on a mocked test. It moves on these numbers.

```bash
# A handful of figures with known panel structure. The gold corpuscle is
# ideal; any built corpuscle with multi-panel plates will do.
cd <corpuscle>

for DT in float32 bfloat16 float16; do
  echo "=== $DT ==="
  # Peak RSS and wall clock, plus the ROIs actually produced.
  CORPUS_VLM_DTYPE=$DT /usr/bin/time -l \
    corpus run --only vision --figure-panels vision-local \
    2>&1 | tee "vision-$DT.log"
  # The loader prints the dtype and its estimated footprint up front:
  #   local VLM dtype=bfloat16 (bfloat16), ~15.2 GB of weights ...
  grep -E "local VLM dtype|maximum resident" "vision-$DT.log"
done
```

Then compare, in this order — the second question matters more than the first:

1. **Does it load and run at all?** float32 may not fit; note the machine's RAM.
2. **Are the ROIs the same?** This is the one that decides it. Compare
   `rois` across the three runs for the same figures — count, and boxes to
   within a few pixels. Half precision that quietly finds fewer or sloppier
   panels is worse than float32 that does not fit, because the corpuscle
   still builds and nothing flags it.
3. **Peak memory and wall clock**, for the record.

Report all three in [#258](https://github.com/caseywdunn/corpus/issues/258).
If bf16 and fp16 both match float32's ROIs, the default becomes half
precision on MPS and the knob stays for the exception. If they differ, the
comment the issue asks for gets written instead — naming the op that fails,
so the next person does not re-litigate it.

**Worth the same sitting:** the `docling==2.94.0` pin
([#98](https://github.com/caseywdunn/corpus/issues/98) follow-up) has been
waiting on the same hardware. v1.2's fidelity harness now gives it a criterion
it never had — score 2.95/2.96 against the gold set rather than against
impressions.

## (2) linux-x86_64 — Bouchet (clean env)

```bash
module load miniconda                      # YCRC convention

# Clean env recreate — the whole point of this leg.
conda env remove -n corpus --yes
conda env create -f environment.yaml
conda activate corpus
pip install -e .
bash tools/install_tessdata.sh

# Grobid via Apptainer/Singularity — see INSTALL.md#grobid-image-choice-and-hosts-without-docker
# On the cluster proper, submit slurm/batch_grobid.sh instead (dev_docs/BOUCHET.md
# §6); it does this correctly and keeps the service off the login node. The two
# writable binds below are not optional — without the tmp bind Grobid answers
# HTTP 500 to every request, and without the logs bind it crashes on startup.
singularity build --force grobid.sif docker://lfoppiano/grobid:0.8.1
GROBID_TMP=$(mktemp -d) && GROBID_LOGS=$(mktemp -d)
singularity run --pwd /opt/grobid --bind $HOME \
  --bind "$GROBID_TMP:/opt/grobid/grobid-home/tmp" \
  --bind "$GROBID_LOGS:/opt/grobid/logs" grobid.sif &

# Models take ~30-60 s to load; the backgrounded run above has not bound the
# port yet, so poll instead of probing once.
until [ "$(curl -fs http://localhost:8070/api/isalive 2>/dev/null)" = true ]; do
  sleep 5; echo "waiting for grobid..."
done                                       # expect: true

cd demo && corpus -v check                 # -v required to see the ok lines
corpus -v run --no-vision                  # wall time depends on Bouchet
                                           # load + WoRMS API rate; budget
                                           # 30–45 min for the demo.
corpus status --report                     # expect: 4 / 4 done
jq '.paper_count' output/corpus_bundle/bundle_manifest.json   # expect: 4
python tools/smoke_test_sse.py demo/output/corpus_bundle --port 18080
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
[INSTALL.md](../INSTALL.md#apple-silicon-arm64-native-conda-required).

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

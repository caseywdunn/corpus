"""Grobid must outlive stage 1, and the walls are what guarantee it.

Grobid is submitted *before* stage 1. If it dies first, the papers stage 1
processes after that get placeholder metadata — and implicit resume will not
retry them, because their inputs are unchanged. The loss is silent, which is
the worst shape a data bug can have.

That invariant used to be bought by putting Grobid on `week` (48 h) while
stage 1 ran on `day` (24 h). `week` became unusable — on 2026-09-08 it had 20
idle CPUs of 1792 and a four-day estimated wait for a 4-CPU job, while `day`
had ~2000 idle — so both now run on `day` and the ordering is bought by the
walltimes instead: Grobid 24 h, stage 1 20 h.

Two numbers in two files now carry a silent-data-loss invariant. This test is
the thing that notices when someone raises one of them.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

SLURM = Path(__file__).resolve().parent.parent / "slurm"
GROBID = SLURM / "batch_grobid.sh"
STAGE1 = SLURM / "batch_process_corpus.sh"


def _directive(script: Path, name: str) -> str:
    pattern = re.compile(rf"^#SBATCH\s+--{name}=(\S+)", re.MULTILINE)
    found = pattern.findall(script.read_text())
    assert found, f"no --{name} in {script.name}"
    # Last wins, the way sbatch reads them.
    return found[-1]


def _seconds(walltime: str) -> int:
    """SLURM walltime to seconds. Handles `D-HH:MM:SS` and `HH:MM:SS`."""
    days, _, rest = walltime.rpartition("-")
    parts = [int(p) for p in rest.split(":")]
    while len(parts) < 3:
        parts.append(0)
    h, m, s = parts
    return int(days or 0) * 86400 + h * 3600 + m * 60 + s


# ── the invariant ──────────────────────────────────────────────────────

def test_grobid_outlives_stage_one():
    grobid = _seconds(_directive(GROBID, "time"))
    stage1 = _seconds(_directive(STAGE1, "time"))
    assert grobid > stage1, (
        f"Grobid asks for {grobid}s and stage 1 for {stage1}s. Grobid is "
        f"submitted first, so an equal-or-shorter wall means it dies first "
        f"and the tail of the corpus silently gets placeholder metadata that "
        f"resume will never retry. Raise Grobid's wall or lower stage 1's."
    )


def test_the_gap_is_wide_enough_to_be_deliberate():
    """A one-minute margin would satisfy the ordering and mean nothing.
    Grobid also has to survive its own startup and model load."""
    gap = _seconds(_directive(GROBID, "time")) - _seconds(_directive(STAGE1, "time"))
    assert gap >= 3600, f"only {gap}s of headroom between Grobid and stage 1"


# ── the walls have to be legal where the jobs are sent ─────────────────

# Bouchet rejects a `week` submission shorter than 24 h, pointing at `day`
# or `scavenge`; `day` caps at exactly 24 h. So a wall is only valid for
# some partitions, and the pair above only works on `day`.
_PARTITION_MAX_SECONDS = {"day": 24 * 3600, "week": 7 * 86400}
_PARTITION_MIN_SECONDS = {"week": 24 * 3600}


@pytest.mark.parametrize("script", [GROBID, STAGE1], ids=lambda p: p.name)
def test_the_walltime_fits_the_partition(script):
    partition = _directive(script, "partition")
    wall = _seconds(_directive(script, "time"))
    ceiling = _PARTITION_MAX_SECONDS.get(partition)
    if ceiling is not None:
        assert wall <= ceiling, (
            f"{script.name} asks {partition} for {wall}s; it caps at {ceiling}s"
        )
    floor = _PARTITION_MIN_SECONDS.get(partition)
    if floor is not None:
        assert wall >= floor, (
            f"{script.name} asks {partition} for {wall}s, below its {floor}s "
            f"minimum — Bouchet rejects such a submission outright"
        )


def test_grobid_is_not_parked_on_a_partition_it_cannot_get():
    """Not a style preference: `week` is where this job spent four days
    waiting behind the consumer that depends on it."""
    assert _directive(GROBID, "partition") != "week"


def test_stage_one_and_grobid_share_a_partition():
    """They no longer trade off against each other's queues. If these
    diverge again, the ordering argument above stops being about walls."""
    assert _directive(GROBID, "partition") == _directive(STAGE1, "partition")

"""SLURM scripts use `module reset`, not `module purge` (#252).

`sbatch --export=ALL` propagates LOADEDMODULES / _LMFILES_, so a batch job
starts believing miniconda is already loaded. `module purge` then unloads it
and the modulefile's unload hook calls `conda` — a shell function that is
not exported and does not exist in a non-interactive shell. Every job
therefore opened its stderr with two alarming lines:

    /apps/.../init/bash: line 169: conda: command not found
    CondaError: Run 'conda init' before 'conda deactivate'

The jobs ran fine. The noise did not: during a failed Stage 1 launch the
real cause was a missing taxonomy.sqlite (#251), and `conda: command not
found` sat above it and drew the first hypothesis.

`module purge` also does not do what its presence implies — StdEnv is sticky
and survives it. `module reset` restores that same sticky default, matches
YCRC's documented convention, and emits one informational line.
"""
from __future__ import annotations

from pathlib import Path

import pytest

SLURM = Path(__file__).resolve().parent.parent / "slurm"
SCRIPTS = sorted(SLURM.glob("batch_*.sh"))


def test_there_are_scripts_to_check():
    assert SCRIPTS, "no slurm/batch_*.sh found; this guard would pass vacuously"


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_no_module_purge(script: Path):
    for i, line in enumerate(script.read_text().splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue          # the explanatory comment names it on purpose
        assert "module purge" not in stripped, (
            f"{script.name}:{i} uses `module purge`; use `module reset` (#252)"
        )

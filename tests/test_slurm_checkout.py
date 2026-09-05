"""Build jobs must load the selected checkout despite an older installation."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


REPO = Path(__file__).resolve().parent.parent
PATHS = REPO / "slurm" / "bouchet_paths.sh"


@pytest.fixture
def checkout_env(tmp_path):
    stale = tmp_path / "old-install"
    for name in ("pipeline", "bib", "mcpsrv"):
        package = stale / name
        package.mkdir(parents=True)
        (package / "__init__.py").write_text("")
    env = {
        **os.environ,
        "REPO_DIR": str(REPO),
        "CORPUS_CONFIG": str(tmp_path / "config.yaml"),
        "PYTHONPATH": str(stale),
        "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"],
    }
    env.pop("CORPUS_BUILD_GIT_SHA", None)
    return env


def test_selected_checkout_reaches_console_script_and_phase_subprocess(
    tmp_path, checkout_env,
):
    # Like an installed `corpus` script: sys.path[0] is its bin directory.
    entrypoint = tmp_path / "bin" / "corpus-probe.py"
    entrypoint.parent.mkdir()
    probe = (
        "import importlib.util, json; "
        "print(json.dumps({n: importlib.util.find_spec(n).origin "
        "for n in ('pipeline', 'bib', 'mcpsrv')}))"
    )
    entrypoint.write_text(
        "import os, subprocess, sys\n"
        + probe + "\n"
        + f"subprocess.run([sys.executable, '-c', {probe!r}], "
        "cwd=os.path.join(os.environ['REPO_DIR'], 'pipeline'), check=True)\n"
    )
    before = subprocess.run(
        [sys.executable, str(entrypoint)], env=checkout_env,
        cwd=tmp_path, text=True, capture_output=True, check=True,
    )
    assert all("old-install" in line for line in before.stdout.splitlines())

    after = subprocess.run(
        ["bash", "-eu", "-c",
         'source "$1"; corpus_check_checkout; "$2" "$3"',
         "probe", str(PATHS), sys.executable, str(entrypoint)],
        env=checkout_env, cwd=tmp_path, text=True, capture_output=True, check=True,
    )
    records = [json.loads(line) for line in after.stdout.splitlines()
               if line.startswith("{")]
    assert len(records) == 2
    for record in records:
        assert record == {
            name: str(REPO / name / "__init__.py")
            for name in ("pipeline", "bib", "mcpsrv")
        }


def test_preflight_rejects_import_path_changed_after_activation(tmp_path, checkout_env):
    result = subprocess.run(
        ["bash", "-eu", "-c",
         'source "$1"; export PYTHONPATH="$2"; corpus_check_checkout',
         "probe", str(PATHS), checkout_env["PYTHONPATH"]],
        env=checkout_env, cwd=tmp_path, text=True, capture_output=True,
    )
    assert result.returncode != 0
    assert "Checkout mismatch" in result.stderr


def test_preflight_rejects_commit_changed_since_submission(tmp_path, checkout_env):
    result = subprocess.run(
        ["bash", "-eu", "-c", 'source "$1"; corpus_check_checkout',
         "probe", str(PATHS)],
        env={**checkout_env, "CORPUS_BUILD_GIT_SHA": "0" * 40},
        cwd=tmp_path, text=True, capture_output=True,
    )
    assert result.returncode != 0
    assert "Checkout changed since submission" in result.stderr

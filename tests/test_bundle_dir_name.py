"""The distilled bundle is `corpus_bundle/`, and `_serve/` still works (#273).

`_serve` was the one directory in a corpuscle designed to be moved away
from the build that produced it, and the one whose name said nothing
about what it is. Landed in an S3 bucket or beside three sibling bundles
it identified neither the project nor the artifact — and the leading
underscore said the opposite of the truth, since by convention `_foo`
reads as private scratch you may delete and this is the only deliverable
in the tree.

The rename cannot be a flag day. Existing corpuscles have a `_serve/`,
and MCP clients, operator scripts and S3 prefixes are pointed at it. So
an existing one is read *and updated in place*, because a rename that
only ever wrote the new name would leave a stale bundle behind for a
client to keep serving — which is worse than an ugly directory name.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from mcpsrv.bundle import (
    BUNDLE_DIR_NAME,
    LEGACY_BUNDLE_DIR_NAME,
    is_bundle_dir,
    resolve_bundle_dir,
)


def _manifest(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    (path / "bundle_manifest.json").write_text(json.dumps({"bundle_version": "1"}))
    return path


# --- the names ----------------------------------------------------------


def test_the_new_name_is_self_describing_and_has_no_underscore():
    assert BUNDLE_DIR_NAME == "corpus_bundle"
    assert not BUNDLE_DIR_NAME.startswith("_")
    assert LEGACY_BUNDLE_DIR_NAME == "_serve"


# --- resolution ---------------------------------------------------------


def test_a_fresh_corpuscle_gets_the_new_name(tmp_path):
    path, is_legacy = resolve_bundle_dir(tmp_path)
    assert path == tmp_path / BUNDLE_DIR_NAME
    assert is_legacy is False


def test_an_existing_legacy_bundle_is_used_and_flagged(tmp_path):
    _manifest(tmp_path / LEGACY_BUNDLE_DIR_NAME)
    path, is_legacy = resolve_bundle_dir(tmp_path)
    assert path == tmp_path / LEGACY_BUNDLE_DIR_NAME
    assert is_legacy is True


def test_the_new_name_wins_when_both_exist(tmp_path):
    """Not an expected layout, but if a corpuscle has both, the one the
    current code writes is the current one."""
    _manifest(tmp_path / LEGACY_BUNDLE_DIR_NAME)
    _manifest(tmp_path / BUNDLE_DIR_NAME)
    path, is_legacy = resolve_bundle_dir(tmp_path)
    assert path == tmp_path / BUNDLE_DIR_NAME
    assert is_legacy is False


def test_the_demo_fixture_still_resolves(tmp_path):
    """The checked-in demo corpuscle predates the rename; it must keep
    working without being rebuilt."""
    root = tmp_path / "output"
    _manifest(root / LEGACY_BUNDLE_DIR_NAME)
    assert resolve_bundle_dir(root)[0].name == LEGACY_BUNDLE_DIR_NAME


# --- identifying a bundle -----------------------------------------------


def test_a_bundle_is_identified_by_its_manifest_not_its_name(tmp_path):
    """#273's own argument: the basename check was already redundant
    beside the manifest check and should not be replaced by a check
    against the new name. A bundle that has been renamed or relocated is
    still a bundle — which is the point of a portable name."""
    for name in ("corpus_bundle", "_serve", "siphonophores-v1.4", "bundle"):
        assert is_bundle_dir(_manifest(tmp_path / name)), name


def test_a_directory_without_a_manifest_is_not_a_bundle(tmp_path):
    (tmp_path / BUNDLE_DIR_NAME).mkdir()
    assert not is_bundle_dir(tmp_path / BUNDLE_DIR_NAME)
    assert not is_bundle_dir(tmp_path / "nonexistent")


def test_the_basename_check_is_gone_from_serve_check():
    """It degraded safely, but it was dead weight that would silently
    stop contributing on any rename."""
    import inspect

    # `from mcpsrv import main` resolves the lazy accessor in
    # mcpsrv/__init__.py to the *function*, not the module.
    from mcpsrv.main import _serve_check
    src = inspect.getsource(_serve_check)
    assert 'build_dir.name == "_serve"' not in src
    assert f'build_dir.name == "{BUNDLE_DIR_NAME}"' not in src, (
        "the basename check must not be reintroduced under the new name"
    )
    assert "is_bundle_dir(build_dir)" in src


# --- what writes it -----------------------------------------------------


def test_the_distiller_updates_a_legacy_bundle_in_place(tmp_path, monkeypatch, capsys):
    """Two bundles in one corpuscle is how a client pointed at the old
    path ends up serving a stale one."""
    from pipeline import cli

    _manifest(tmp_path / LEGACY_BUNDLE_DIR_NAME)
    seen = {}

    def _run(cmd, **_kw):
        seen["cmd"] = cmd
        return type("R", (), {"returncode": 0})()

    monkeypatch.setattr(cli.subprocess, "run", _run)
    cli._distill_bundle(tmp_path)
    assert str(tmp_path / LEGACY_BUNDLE_DIR_NAME) in seen["cmd"]
    out = capsys.readouterr().out
    assert "pre-1.4 name" in out
    assert "mv " in out, "the migration command should be in the nudge"


def test_the_distiller_writes_the_new_name_for_a_fresh_corpuscle(
    tmp_path, monkeypatch, capsys,
):
    from pipeline import cli

    seen = {}

    def _run(cmd, **_kw):
        seen["cmd"] = cmd
        return type("R", (), {"returncode": 0})()

    monkeypatch.setattr(cli.subprocess, "run", _run)
    cli._distill_bundle(tmp_path)
    assert str(tmp_path / BUNDLE_DIR_NAME) in seen["cmd"]
    assert "pre-1.4 name" not in capsys.readouterr().out


# --- nothing still hardcodes the old name where it matters --------------


def test_no_shipped_script_or_workflow_hardcodes_the_old_path():
    """A fresh build writes the new name, so a leftover `_serve` path in
    CI, deploy or the SLURM scripts is a break, not a compatibility
    nicety."""
    import re

    root = Path(__file__).resolve().parent.parent
    pattern = re.compile(r"(?<![A-Za-z0-9_])_serve(?![A-Za-z0-9_])")
    offenders = []
    for sub in (".github/workflows", "deploy", "slurm", "tools"):
        for path in (root / sub).rglob("*"):
            if not path.is_file() or path.suffix not in (".yml", ".sh", ".py", ".conf"):
                continue
            for i, line in enumerate(path.read_text().splitlines(), 1):
                if pattern.search(line):
                    offenders.append(f"{path.relative_to(root)}:{i}: {line.strip()}")
    assert not offenders, "still pointing at the old bundle name:\n" + "\n".join(offenders)

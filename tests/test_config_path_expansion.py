"""Config paths expand `~` and `$VARS`, and an unset variable stays literal.

YAML has no shell, so every expansion a user expects has to be done here. `~`
was already handled; `$VARS` matters because it is what lets one config travel
between machines — `output_dir: $CORPUS_DATA/corpuscles/viburnum_20260913`
resolves to a workstation's data root and a cluster's without editing the file.

That portability is not cosmetic. A corpuscle is derived data several times the
size of its sources, so it belongs wherever a given host keeps large things,
which is rarely the same path twice and rarely inside the library's git repo.
Without expansion, the only ways to express that are an absolute path that is
wrong on every other machine, or a relative one that buries gigabytes in a
checkout.

The unset case is the one worth pinning. Shell would turn `$NOPE/corpuscles/x`
into `/corpuscles/x` — a silent write to the filesystem root, as root-adjacent a
failure as this code can produce. `os.path.expandvars` leaves it literal
instead, so it fails as a path containing `$NOPE`, which names the actual
mistake. This test exists so nobody "fixes" that into shell semantics.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from pipeline.cli import _resolve_against


@pytest.fixture
def config_path(tmp_path: Path) -> Path:
    cfg = tmp_path / "corpuscle" / "config.yaml"
    cfg.parent.mkdir(parents=True)
    cfg.write_text("", encoding="utf-8")
    return cfg


def test_environment_variable_is_expanded(config_path, monkeypatch):
    monkeypatch.setenv("CORPUS_DATA", "/data/corpus")
    result = _resolve_against(config_path, Path("$CORPUS_DATA/corpuscles/x"))
    assert result == Path("/data/corpus/corpuscles/x")


def test_braced_form_is_expanded(config_path, monkeypatch):
    """`${VAR}` is the form a reader writes next to other path segments."""
    monkeypatch.setenv("CORPUS_DATA", "/data/corpus")
    result = _resolve_against(config_path, Path("${CORPUS_DATA}/corpuscles/x"))
    assert result == Path("/data/corpus/corpuscles/x")


def test_unset_variable_is_left_literal_not_emptied(config_path, monkeypatch):
    """The safety property: never silently resolve to the filesystem root."""
    monkeypatch.delenv("CORPUS_NOT_SET", raising=False)
    result = _resolve_against(config_path, Path("$CORPUS_NOT_SET/corpuscles/x"))

    assert "$CORPUS_NOT_SET" in str(result), (
        "an unset variable must stay literal so the error names it"
    )
    assert result != Path("/corpuscles/x"), (
        "shell semantics here would write to the filesystem root"
    )


def test_tilde_still_expands(config_path):
    result = _resolve_against(config_path, Path("~/data/pdfs"))
    assert result == Path.home() / "data" / "pdfs"


def test_relative_paths_still_resolve_against_the_config(config_path):
    """The original contract: `cd demo && corpus run` and `--config demo/...`
    from anywhere must agree."""
    result = _resolve_against(config_path, Path("./library"))
    assert result == (config_path.parent / "library").resolve()


def test_absolute_paths_are_untouched(config_path):
    result = _resolve_against(config_path, Path("/srv/corpora/thing"))
    assert result == Path("/srv/corpora/thing")


def test_none_stays_none(config_path):
    assert _resolve_against(config_path, None) is None


def test_expanded_variable_is_treated_as_absolute(config_path, monkeypatch):
    """An expanded value that is absolute must not then be glued onto the
    corpuscle root — that was the bug the `~` handling exists to prevent, and
    it would recur identically for `$VARS`."""
    monkeypatch.setenv("CORPUS_DATA", "/data/corpus")
    result = _resolve_against(config_path, Path("$CORPUS_DATA/out"))
    assert not str(result).startswith(str(config_path.parent))

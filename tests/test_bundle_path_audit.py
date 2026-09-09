"""The absolute-path audit does not flag extracted content (#183).

`_audit_no_absolute_paths` flagged PDF glyph names and OCR garbage as
filesystem paths and raised. On a 20,137-paper corpus that killed
`corpus run --only bundle` after ~3 h and ~370,000 files copied, with no
`bundle_manifest.json` written — leaving `_serve/` complete but
unservable. Thirteen junk strings in three documents blocked the whole
corpuscle, and all thirteen were in content fields: 9 in
`chunks[].headings[0]`, 1 in `chunks[].text`, 1 in `references[].title`.

No shape-based rule can fix that. `/Peswme/` and `/scratch/` are the same
shape, and one flagged value was `/Summary/` — a correct section heading
that OCR wrapped in slashes. So content fields are exempt from the shape
rule and checked only against the build's own root, which can be matched
exactly.

The shape rule stays for every other field, because that is the release
gate: it is what catches a path leak in a field no scrubber knows about
yet.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from mcpsrv.bundle import (
    _ABS_PATH_RE,
    _CONTENT_KEYS,
    _audit_no_absolute_paths,
    _walk_keyed_strings,
)


def _bundle(tmp_path: Path, name: str, **files) -> Path:
    root = tmp_path / name
    d = root / "documents" / "aaaaaaaaaaaa"
    d.mkdir(parents=True)
    for fname, payload in files.items():
        (d / f"{fname}.json").write_text(json.dumps(payload))
    return root


# --- the reported failure ------------------------------------------------


@pytest.mark.parametrize("value", [
    "/CP/D2/CP/D0/DD/D7/CT/D7/BA",                       # glyph names
    "/BT/CR/CZ/D2/D3/DB/D0/CT/CS/CV/D1/CT/D2/D8/D7",
    "/Peswme/",                                          # OCR noise
    "/ifl/fl^Ma5Mm^«Mm«aCram",
    "/Summary/",                                         # a real heading
])
def test_content_that_looks_like_a_path_does_not_fail_the_bundle(tmp_path, value):
    """Every one of these still matches the shape rule — that is why the
    exemption is by field, not by pattern."""
    assert _ABS_PATH_RE.match(value), "fixture no longer reproduces the case"
    root = _bundle(tmp_path, value[:4].replace("/", "_"), chunks={
        "chunks": [{"chunk_id": "c0", "headings": [value], "text": value}],
    })
    assert _audit_no_absolute_paths(root, build_roots=[tmp_path]) == []


def test_a_reference_title_is_content_too(tmp_path):
    root = _bundle(tmp_path, "refs", references={
        "references": [{"title": "/ifl/fl^Ma5Mm^«Mm«aCram"}],
    })
    assert _audit_no_absolute_paths(root, build_roots=[tmp_path]) == []


# --- what must still fail ------------------------------------------------


def test_a_real_leak_in_a_path_field_still_fails(tmp_path):
    root = _bundle(tmp_path, "leak", figures={
        "figures": [{"figure_id": "f1",
                     "file_path": "/nfs/roberts/project/build/f1.png"}],
    })
    offenders = _audit_no_absolute_paths(root, build_roots=[tmp_path])
    assert len(offenders) == 1
    where, value = offenders[0]
    assert where.endswith("figures.json: figures[0].file_path")
    assert value.startswith("/nfs/roberts")


def test_a_path_field_no_scrubber_knows_about_still_fails(tmp_path):
    """The gate's whole purpose. A denylist of content keys means a field
    added later defaults to being checked; an allowlist of path keys would
    have traded this loud failure for a silent miss."""
    root = _bundle(tmp_path, "new", summary={
        "some_field_added_next_year": "/scratch/user/out/",
    })
    assert len(_audit_no_absolute_paths(root, build_roots=[tmp_path])) == 1


def test_the_report_names_the_field_not_just_the_file(tmp_path):
    """The old message told the operator to fix `_scrub_summary` /
    `_scrub_figures` for values in `chunks[].headings`, sending them into
    the wrong subsystem entirely."""
    root = _bundle(tmp_path, "ptr", summary={
        "steps": [{"stage": "x"}, {"stage": "y", "out": "/scratch/z/"}],
    })
    where, _ = _audit_no_absolute_paths(root, build_roots=[tmp_path])[0]
    assert where.endswith("summary.json: steps[1].out")


# --- the build-root rule for content -------------------------------------


def test_a_build_path_inside_content_warns_but_does_not_fail(tmp_path, caplog):
    """Extraction should not inject one, so it is worth reporting — but it
    cannot be scrubbed, and failing a bundle after the copy is the harm
    this issue is about."""
    root = _bundle(tmp_path, "inject", chunks={
        "chunks": [{"chunk_id": "c0", "text": f"see {tmp_path}/build/x here"}],
    })
    with caplog.at_level("WARNING"):
        assert _audit_no_absolute_paths(root, build_roots=[tmp_path]) == []
    assert "not scrubbable, not fatal" in caplog.text
    assert "chunks[0].text" in caplog.text


def test_no_build_roots_means_no_content_check(tmp_path):
    """Callers that do not know the build root still get the shape gate."""
    root = _bundle(tmp_path, "noroot", chunks={
        "chunks": [{"chunk_id": "c0", "text": "/Peswme/"}],
    })
    assert _audit_no_absolute_paths(root) == []


# --- the walker ---------------------------------------------------------


def test_the_keyed_walk_reports_a_usable_pointer():
    tree = {"chunks": [{"chunk_id": "c0", "headings": ["A", "B"]}]}
    found = dict((p, v) for p, _k, v in _walk_keyed_strings(tree))
    assert found["chunks[0].chunk_id"] == "c0"
    assert found["chunks[0].headings[0]"] == "A"
    assert found["chunks[0].headings[1]"] == "B"


def test_list_elements_inherit_their_containing_key():
    """`headings` is a list of strings, and it was the biggest offender —
    9 of the 13. Its elements have to be recognised as content."""
    _p, key, _v = next(iter(_walk_keyed_strings({"headings": ["/Summary/"]})))
    assert key == "headings"
    assert key in _CONTENT_KEYS


def test_the_path_bearing_keys_are_not_exempt():
    """A regression guard on the denylist itself."""
    for key in ("file_path", "filename", "files_created", "path",
                "output_dir", "input_dir"):
        assert key not in _CONTENT_KEYS, key

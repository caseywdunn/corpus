"""CJK whitespace is normalized for comparison, not for fingerprinting (#280).

Two solo builds of the same 35-document gold set — same commit, same
machine, no contention, no OCR timeouts, identical quality flags —
differed on exactly two documents, both Japanese scans, and the
differences were whitespace segmentation only: same glyphs, different
space placement, 186 lines of it.

Japanese does not delimit words with spaces, so where OCR puts them
between CJK characters carries no information. But it changes the digest,
and it buried whatever real difference the acceptance run was looking
for.

**The `--jobs` hypothesis is not supported.** The issue named `jobs=12`
as the likely mechanism. Tested directly on both affected documents —
Yamamori2014 (jpn+eng) and Kawamura1911a pages 14-25 (jpn_vert) — the
OCR output is byte-identical across repeated runs at `--jobs` 1, 4 and
12, and across `OMP_THREAD_LIMIT` 1, 4, 12 and unset. docling is
likewise deterministic on identical input across three runs. So pinning
`--jobs 1` for comparison builds would trade real build time for
nothing, and that remedy should not be re-attempted on this reasoning.
"""
from __future__ import annotations

import pytest

from tools.qc.build_reference import _logical, normalize_cjk_spacing


# --- the observed diffs -------------------------------------------------


@pytest.mark.parametrize("a,b", [
    ("和歌 山県田辺 湾 の海岸 線付近", "和歌 山県 田辺 湾 の 海岸 線 付近"),
    ("田辺 次は ,", "田辺 次 は ,"),
    ("海岸 か ら 表 層", "海岸 から 表 層"),
])
def test_the_real_diffs_collapse_to_the_same_string(a, b):
    """Verbatim from the acceptance run's diff output."""
    assert normalize_cjk_spacing(a) == normalize_cjk_spacing(b)


def test_whitespace_between_cjk_characters_is_dropped():
    assert normalize_cjk_spacing("海岸 線") == "海岸線"
    assert normalize_cjk_spacing("海岸\n線") == "海岸線"
    assert normalize_cjk_spacing("海岸  \t 線") == "海岸線"


# --- what must stay visible ---------------------------------------------


def test_a_latin_word_boundary_is_a_real_difference():
    """Scoped to between two CJK characters precisely so the comparison
    stays sharp everywhere else — a lost space in Latin text is content
    loss, not segmentation noise."""
    assert normalize_cjk_spacing("one two") != normalize_cjk_spacing("onetwo")
    assert normalize_cjk_spacing("one two") == "one two"


def test_a_cjk_to_latin_boundary_keeps_its_space():
    assert normalize_cjk_spacing("海岸 line ここ") == "海岸 line ここ"


def test_a_space_before_punctuation_is_kept():
    """`田辺 次 は ,` keeps its space before the comma, which is not CJK
    — so a punctuation-spacing change would still show up."""
    assert normalize_cjk_spacing("田辺 次 は ,") == "田辺次は ,"


def test_a_dropped_cjk_character_is_still_a_difference():
    """Normalization must not hide lost content."""
    assert normalize_cjk_spacing("和歌 山県") != normalize_cjk_spacing("和歌 山")


@pytest.mark.parametrize("value", ["", None])
def test_empty_and_none_are_passed_through(value):
    assert normalize_cjk_spacing(value) == value


def test_kana_and_hangul_count_as_cjk():
    assert normalize_cjk_spacing("ここ から") == "ここから"      # hiragana
    assert normalize_cjk_spacing("カタ カナ") == "カタカナ"      # katakana
    assert normalize_cjk_spacing("한국 어") == "한국어"          # hangul


# --- where it is applied, and where it is not ---------------------------


def test_it_is_applied_through_the_logical_snapshot(tmp_path):
    """The harness digests `_logical`, so normalization has to sit there
    for the two builds to compare equal."""
    a = _logical({"text": "和歌 山県田辺 湾", "pages": 4}, tmp_path)
    b = _logical({"text": "和歌 山県 田辺 湾", "pages": 4}, tmp_path)
    assert a == b


def test_it_reaches_nested_strings(tmp_path):
    a = _logical({"chunks": [{"text": "海岸 線"}]}, tmp_path)
    b = _logical({"chunks": [{"text": "海岸線"}]}, tmp_path)
    assert a == b


def test_it_does_not_touch_the_pipeline_fingerprint():
    """A fingerprint decides what re-runs. Rewriting the text it hashes
    would change resume behaviour to buy a property only the acceptance
    harness needs — the issue offered that; this is the free half."""
    import inspect

    from pipeline import stages
    src = inspect.getsource(stages)
    assert "normalize_cjk_spacing" not in src


def test_the_normalizer_lives_in_the_harness_not_the_pipeline():
    """`tools/` may read the product; the product may not depend on its
    own measurement tooling (AGENTS.md, enforced by
    tests/test_import_direction.py)."""
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent
    for sub in ("pipeline", "mcpsrv", "bib"):
        for path in (root / sub).rglob("*.py"):
            assert "normalize_cjk_spacing" not in path.read_text(), path

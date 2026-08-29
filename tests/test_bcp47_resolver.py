"""The public BCP-47 → Tesseract pack resolver (#215).

This table's job is to be the *only* copy. A library's annotation pass
resolves its per-document language into an `ocrlang` bib directive and needs
exactly the mapping `scan.py` uses; a local copy in the library repo is how
the two come to disagree with nothing noticing. So these tests pin the
contract that an out-of-repo consumer relies on, including the parts that
would look like oversights to someone tidying up later.
"""
from __future__ import annotations

import pytest

from pipeline.scan import _parse_bcp47, bcp47_to_tesseract


@pytest.mark.parametrize(
    "tag,expected",
    [
        # Bare ISO 639-1 is already a well-formed BCP-47 tag, so everything
        # langdetect emits keeps resolving unchanged.
        ("ru", ["rus"]),
        ("de", ["deu"]),
        ("en", ["eng"]),
        ("la", ["lat"]),
        # The cases ISO cannot name, which are the reason for the issue.
        ("de-Latf", ["deu_latf", "deu"]),
        ("zh-Hans", ["chi_sim"]),
        ("zh-Hant", ["chi_tra"]),
        ("sr-Latn", ["srp_latn"]),
        ("sr-Cyrl", ["srp"]),
        # 639-3 for a language with no 639-1 code.
        ("grc", ["grc"]),
        # Unknown, and the "no linguistic content" tags, resolve to nothing
        # rather than guessing.
        ("xx", []),
        ("und", []),
        ("zxx", []),
        ("", []),
    ],
)
def test_resolves(tag, expected):
    assert bcp47_to_tesseract(tag) == expected


def test_langdetect_legacy_region_tags_still_resolve():
    """`zh-cn` / `zh-tw` encode script as region and are not really ISO.

    They are still what `scan_detection.json` records for Chinese, so a
    consumer comparing against a detection record must not fall off the map.
    """
    assert bcp47_to_tesseract("zh-cn") == ["chi_sim"]
    assert bcp47_to_tesseract("zh-tw") == ["chi_tra"]


def test_subtag_case_is_not_significant():
    """BCP-47 conventions are conventions; tags arrive hand-written."""
    assert bcp47_to_tesseract("de-latf") == bcp47_to_tesseract("de-Latf")
    assert bcp47_to_tesseract("DE-LATF") == bcp47_to_tesseract("de-Latf")
    assert bcp47_to_tesseract("zh-hant") == bcp47_to_tesseract("zh-Hant")
    assert bcp47_to_tesseract("de_Latf") == bcp47_to_tesseract("de-Latf")


def test_script_beats_region_and_trailing_subtags_are_ignored():
    """`de-Latf-DE` is Fraktur; the region and any variant say nothing."""
    assert bcp47_to_tesseract("de-Latf-DE") == ["deu_latf", "deu"]
    assert bcp47_to_tesseract("de-DE-1901") == ["deu"]
    assert bcp47_to_tesseract("zh-Hant-HK") == ["chi_tra"]


def test_parse_identifies_subtags_by_shape():
    assert _parse_bcp47("de-Latf-DE") == ("de", "Latf", "DE")
    assert _parse_bcp47("es-419") == ("es", None, "419")
    assert _parse_bcp47("en") == ("en", None, None)


def test_bare_zh_returns_both_han_packs():
    """A bare `zh` names two scripts and no single pack.

    Both is the honest answer, and it matches what the OSD branch already
    does for a `Han` verdict. Anything narrower would be a coin flip
    presented as a decision.
    """
    assert bcp47_to_tesseract("zh") == ["chi_sim", "chi_tra"]


def test_no_vertical_companions_are_ever_returned():
    """Vertical layout is typesetting, not language or script (#196).

    `jpn_vert` is also measurably *not* safe to union with `jpn`: on the gold
    set's vertical pages `jpn_vert` alone scores 0.574, `jpn` alone 0.246, and
    the union 0.186 — worse than either. Someone making the CJK entries
    symmetric with the Fraktur one would regress that, so assert it directly.
    """
    for tag in ("ja", "ja-JP", "zh", "zh-Hant", "ko", "ko-KR"):
        assert not any(p.endswith("_vert") for p in bcp47_to_tesseract(tag)), tag


def test_result_is_a_fresh_list_the_caller_may_mutate():
    first = bcp47_to_tesseract("de-Latf")
    first.append("frk")
    assert bcp47_to_tesseract("de-Latf") == ["deu_latf", "deu"]


def test_every_pack_named_is_one_the_pipeline_knows():
    """No entry may name a pack the rest of `scan.py` has never heard of.

    The failure this catches is a typo in a table nothing else reads —
    written into a bib, then silently dropped at availability time on every
    host, looking like a missing install rather than a bad mapping.
    """
    from pipeline.scan import (
        _BCP47_AMBIGUOUS,
        _BCP47_LANG_PACKS,
        _BCP47_REGION_PACKS,
        _BCP47_SCRIPT_PACKS,
        _ISO_TO_TESSERACT,
        _SCRIPT_PACK_FAMILIES,
        _SCRIPT_TO_TESSERACT,
        _VERTICAL_COMPANION,
    )

    known = set(_ISO_TO_TESSERACT.values())
    known |= {p for packs in _SCRIPT_PACK_FAMILIES.values() for p in packs}
    known |= {p for packs in _SCRIPT_TO_TESSERACT.values() for p in packs}
    known |= set(_VERTICAL_COMPANION.values())

    named = set()
    for table in (_BCP47_SCRIPT_PACKS, _BCP47_REGION_PACKS,
                  _BCP47_LANG_PACKS, _BCP47_AMBIGUOUS):
        for packs in table.values():
            named.update(packs)

    assert named <= known, f"unknown pack names: {sorted(named - known)}"

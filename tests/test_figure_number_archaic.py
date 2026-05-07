"""Tests for archaic-format figure-number parsing (#16).

The 19th-c./early-20th-c. tail of the corpus uses caption conventions
parse_figure_number didn't originally cover:

- German "Tafel" (Taf. III.)
- Latin "Tabula"  (Tab. XII.)
- Roman numerals as the figure number itself (Plate IV.)

This test file locks in the new behaviour.
"""
from __future__ import annotations

import pytest

from pipeline.figures import _roman_to_int, parse_figure_number


# ---------------------------------------------------------------------------
# _roman_to_int
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("roman, expected", [
    ("I", 1),
    ("IV", 4),
    ("V", 5),
    ("IX", 9),
    ("X", 10),
    ("XL", 40),
    ("L", 50),
    ("XC", 90),
    ("C", 100),
    ("CD", 400),
    ("D", 500),
    ("CM", 900),
    ("M", 1000),
    ("XII", 12),
    ("XLII", 42),
    ("MCMLXXXIV", 1984),
    ("iv", 4),
    ("xii", 12),
    ("Iv", 4),  # mixed case
])
def test_roman_to_int_valid(roman, expected):
    assert _roman_to_int(roman) == expected


@pytest.mark.parametrize("not_roman", [
    "",
    "5",       # digits aren't Roman
    "XYZ",     # invalid letters
    "IIII",    # invalid form (should be IV)
    "VV",      # invalid (should be X)
    "ABC",     # not Roman
    "10",
])
def test_roman_to_int_rejects_non_romans(not_roman):
    assert _roman_to_int(not_roman) is None


# ---------------------------------------------------------------------------
# parse_figure_number — pre-existing behavior preserved
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("caption, expected", [
    ("Figure 3. Nectophore of Nanomia.", "3"),
    ("Fig. 4. detail.", "4"),
    ("Fig. 4a. multi-panel.", "4a"),
    ("Plate 2.", "2"),
    ("Abb. 5. Querschnitt.", "5"),
    ("Tav. 7. Esemplare.", "7"),
    ("Рис. 8. Нектофор.", "8"),
    ("", None),
    ("Just a paragraph of text.", None),
    ("See Figure 3 in Smith 1995.", None),  # not anchored at start
])
def test_parse_figure_number_legacy(caption, expected):
    assert parse_figure_number(caption) == expected


# ---------------------------------------------------------------------------
# parse_figure_number — archaic formats added in #16
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("caption, expected", [
    # Tafel (German "plate") — common in 19th-c. German monographs
    ("Taf. III. Pelagia.", "3"),
    ("Tafel V. Detail.", "5"),
    ("Taf. 12. modern Arabic.", "12"),
    # Tabula (Latin "plate") — Linnean-era usage
    ("Tab. XII. Animalia.", "12"),
    ("Tabula IV.", "4"),
    # Roman numeral with familiar prefixes
    ("Plate IV. Diphyes.", "4"),
    ("Pl. iv. lowercase.", "4"),
    ("Plate VIII.", "8"),
    ("Figure XL. forty.", "40"),
    ("Fig. III. roman.", "3"),
    # Mixed numbers stay Arabic
    ("Plate 2.", "2"),
    ("Fig. 4.", "4"),
])
def test_parse_figure_number_archaic(caption, expected):
    assert parse_figure_number(caption) == expected


def test_roman_normalized_for_join_with_body_refs():
    """Body text says 'Fig. 4'; caption says 'Plate IV.'. Both
    normalize to '4' so a downstream join key matches.
    """
    assert parse_figure_number("Plate IV. Diphyes.") == "4"
    # body text would be matched by _FIGURE_REF_RE which captures \d+;
    # the join key on both sides is "4".

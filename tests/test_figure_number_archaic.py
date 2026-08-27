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


# ---------------------------------------------------------------------------
# parse_figure_number — OCR-damaged and untranslated openers (#205)
# ---------------------------------------------------------------------------
#
# Measured against 320 stored captions from the 35-document gold corpuscle,
# the damaged spellings of "FIG." are *more common than the correct one*.
# Leading tokens, by count:
#
#     Fic 65   Fig 53   PLATE 35   Figur 17   FIGURE 9   FIG 8
#     Fi   8   FiG   6  Plate  3   Figg   1   Frc    1   Fie 1   Puc 1
#
# Before this, `parse_figure_number` fired for 32.1% of captions; after, 67.6%,
# with precision going *up* from 97.0% to 98.2%.


@pytest.mark.parametrize("caption, expected", [
    # "FIG." set in small caps, misread document-wide. The single largest
    # cause of missing figure numbers in the reference corpus.
    ("Fic.  4 Différenciation du pneumatophore de Nanomia bijuga.", "4"),
    ("Frc.  3 Différenciation du pneumatophore.", "3"),
    ("Fie. 7. Something.", "7"),
    ("Fi. 3 Caption.", "3"),
    # OCR inserting a digit *inside* the opener. Without allowing digits in
    # the opener, the match stops at "Fi" and captures the noise ("1G") as
    # the figure number — which is worse than finding nothing.
    ("Fi1G.  24. Halistemma striata sp.n, from Bermuda.", "24"),
    ("Fi16.  89. Sulculeolaria angusta Totton.", "89"),
    # German "Figur" — not damage at all, simply a spelling the old prefix
    # could not reach, since `fig(?:ure|\.?)` stops before the "ur".
    ("Figur  23. Diphyopsis campanulifera Quoi und Gaimard.", "23"),
    # Italian plural with a range: several figures under one caption. The
    # first number is what joins to a body-text reference.
    ("Figg.  2-5.  -  Sezioni di gonofori di Muggiaca kochi.", "2"),
    # Cyrillic "Рис" read as Latin lookalikes (Р→P, и→u, с→c).
    ("Puc. 7. Нектофор.", "7"),
])
def test_parse_figure_number_ocr_damaged_openers(caption, expected):
    assert parse_figure_number(caption) == expected


@pytest.mark.parametrize("caption, expected", [
    ("Figure 3. Nectophore of Nanomia.", "3"),
    ("FIGURE 5 Caption.", "5"),
    ("Fig. 5. Caption.", "5"),
    ("PLATE XXI", "21"),
    ("Plate IV. Diphyes.", "4"),
    ("Taf. III. Pelagia.", "3"),
    ("Рис. 8. Нектофор.", "8"),
    ("Text-figure 12. Caption.", "12"),
    ("Fig. 3a. Panel.", "3a"),
])
def test_tolerant_opener_does_not_regress_the_correct_spellings(caption, expected):
    """The tolerant opener is an addition, not a replacement — every spelling
    that parsed before must still parse, including Roman normalisation and
    sub-numbering."""
    assert parse_figure_number(caption) == expected


def test_roman_branch_does_not_eat_an_ordinary_word():
    """`[IVXLCDM]` overlaps with real words: "caption" begins with a "c",
    which is Roman 100. The old docstring predicted this; the tolerant opener
    made it reachable, so the Roman branch now requires a non-letter after it.
    """
    assert parse_figure_number("Fig5 caption") == "5"
    assert parse_figure_number("Fig. 5 caption") == "5"
    # A genuine trailing Roman still works when it ends the token.
    assert parse_figure_number("Plate XXI") == "21"
    assert parse_figure_number("Taf. III. Pelagia.") == "3"


@pytest.mark.parametrize("caption", [
    "",
    "Just a paragraph of text.",
    "See Figure 3 in Smith 1995.",          # a reference, not a caption opener
    "[Copied from Totton, 1965, Pl. IX]",   # a citation of another work's plate
    "Differentiation of the pneumatophore.",
])
def test_tolerant_opener_still_requires_a_leading_label(caption):
    """The opener is loose but the anchor is not. A caption must *begin* with
    a label followed by a number; a figure mentioned mid-sentence is a
    cross-reference and binding it would invent a number for the figure."""
    assert parse_figure_number(caption) is None


def test_body_text_reference_matching_was_not_loosened():
    """`_FIGURE_REF_RE` scans running prose, where a fuzzy prefix would bind
    ordinary words to figure numbers. The tolerance is deliberately confined
    to `_FIGURE_NUMBER_IN_CAPTION_RE`, which is anchored to a caption start.
    """
    from pipeline.figures import _FIGURE_REF_RE
    assert _FIGURE_REF_RE.search("as shown in Fig. 4 above")
    # The damaged spellings are NOT accepted in body text.
    assert not _FIGURE_REF_RE.search("as shown in Fic. 4 above")
    assert not _FIGURE_REF_RE.search("we counted for 3 specimens")


# ---------------------------------------------------------------------------
# The tolerant opener must not read ordinary prose as a figure number (#204)
# ---------------------------------------------------------------------------
#
# A regression introduced with the tolerant opener and caught by inspecting a
# promoted image, which turned out to be a handwritten marginal scribble. The
# caption bound to it began "from  the  coasts  of  British  Columbia" — the
# opener matched "fro" and the capture then read the "m" of "from" as Roman
# numeral M, giving it figure number 1000.
#
# The fix separates the two opener classes. A correctly spelled prefix may be
# followed by a Roman numeral with no separator ("PLATE XXI"); an OCR-damaged
# one must be followed by a period. Every damaged spelling in the corpus
# carries the period, so requiring it costs nothing.


@pytest.mark.parametrize("caption", [
    "from  the  coasts  of  British  Columbia «  5 (after  Bigelow)",
    "from 1000 specimens",
    "for the rest of the section",
    "found in shallow water",
])
def test_prose_beginning_with_an_f_word_is_not_a_figure_number(caption):
    assert parse_figure_number(caption) is None


def test_a_damaged_opener_still_needs_its_period():
    """`Fic. 4` is a caption; `Fic 4` without the period is not distinguishable
    from prose, and every damaged spelling observed carries the period."""
    assert parse_figure_number("Fic. 4 Différenciation") == "4"
    assert parse_figure_number("Fic 4 Différenciation") is None


def test_a_correct_prefix_still_needs_no_period():
    """The exact branch keeps its freedom — plate captions are routinely
    written without one."""
    assert parse_figure_number("PLATE XXI") == "21"
    assert parse_figure_number("Figur  23. Diphyopsis") == "23"
    assert parse_figure_number("Figure 5 Caption") == "5"

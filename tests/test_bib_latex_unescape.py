"""BibTeX values are LaTeX source, not plain text (#177).

The reported failure: on a 702-paper build, five titles were served
reading ``... (R.Br.) A.Braun \\& Vatke``. The ``.bib`` was correct —
``\\&`` is the required escaping for a literal ampersand — but corpus
treated the field as plain text on ingestion, so the escape survived into
``metadata.json``, ``_serve/``, and every MCP tool that returns a title.
Any citation a client emitted carried the backslash.

``&`` is only the common case. The same applies to every character a
BibTeX writer must escape (``\\% \\_ \\# \\$``), to group-protection
braces (``{DNA}``, ``{V}iburnum``) which are ubiquitous in real bib files,
and to accent commands (``M{\\"u}ller``, ``Ma\\'nko``) which this corpus is
full of because its authors are largely European.

Both halves are tested: unescaping on the way in, and escaping on the way
out, because ``corpus bib export`` emits a ``.bib`` a user may re-import.
"""
from __future__ import annotations

import pytest

from bib.export import _escape
from bib.parser import bib_entry_to_metadata, parse_bibtex, unescape_bib_value


# --- escaped specials ---------------------------------------------------------

@pytest.mark.parametrize("src,want", [
    (r"Jammu \& Kashmir", "Jammu & Kashmir"),
    (r"(R.Br.) A.Braun \& Vatke", "(R.Br.) A.Braun & Vatke"),
    (r"50\% of taxa", "50% of taxa"),
    (r"gene\_1", "gene_1"),
    (r"\#3", "#3"),
    (r"\$5", "$5"),
])
def test_escaped_specials_are_unescaped(src, want):
    assert unescape_bib_value(src) == want


# --- group-protection braces --------------------------------------------------

@pytest.mark.parametrize("src,want", [
    (r"{DNA} sequencing", "DNA sequencing"),
    (r"{V}iburnum tinus", "Viburnum tinus"),
    (r"The {DNA} of {Nanomia bijuga}", "The DNA of Nanomia bijuga"),
])
def test_group_protection_braces_do_not_reach_metadata(src, want):
    assert unescape_bib_value(src) == want


# --- accent commands ----------------------------------------------------------

@pytest.mark.parametrize("src,want", [
    (r"M{\"u}ller", "Müller"),          # braced diaeresis
    (r'M\"uller', "Müller"),            # bare diaeresis
    (r"Ma\'nko", "Mańko"),              # acute on n — this corpus's own author
    (r"Ma{\'n}ko", "Mańko"),
    (r"Fran\c{c}ois", "François"),      # cedilla, letter command
    (r"\v{S}koda", "Škoda"),            # caron
    (r"Wei\ss", "Weiß"),                # standalone letter command
    (r"\O rsted", "Ørsted"),            # space terminates the command
    (r"\O{}rsted", "Ørsted"),           # empty group terminates it too
])
def test_accents_become_unicode(src, want):
    assert unescape_bib_value(src) == want


def test_letter_command_space_is_a_terminator_not_content():
    """In LaTeX a space *ends* a letter command rather than being text, so
    a .bib that wants the space writes ``\\ss{} and``."""
    assert unescape_bib_value(r"Wei\ss{} and more") == "Weiß and more"
    assert unescape_bib_value(r"Wei\ss and") == "Weißand"


def test_output_is_nfc_normalized():
    """Combining marks are composed, so a served name compares equal to
    one typed directly rather than differing only by normalization form."""
    import unicodedata
    got = unescape_bib_value(r"M{\"u}ller")
    assert got == unicodedata.normalize("NFC", got)
    assert got == "Müller"


# --- nothing gained, nothing lost ---------------------------------------------

@pytest.mark.parametrize("src", [
    "A perfectly ordinary title",
    "Notes on Nanomia bijuga (Delle Chiaje, 1844)",
    "",
])
def test_plain_values_are_untouched(src):
    assert unescape_bib_value(src) == src.strip()


def test_no_stray_backslash_survives_into_metadata():
    """The user-visible symptom in #177: a backslash in a served title."""
    bib = r"""
    @article{braun2020,
      title = {Lineae chinensis (R.Br.) A.Braun \& Vatke and {DNA} barcoding},
      author = {M{\"u}ller, Hans and Ma\'nko, Maciej K.},
      journal = {Journal of 50\% Things},
      year = {2020},
      file = {Braun2020.pdf}
    }
    """
    meta = bib_entry_to_metadata(parse_bibtex(bib)[0], "Braun2020.pdf")
    assert meta["title"] == (
        "Lineae chinensis (R.Br.) A.Braun & Vatke and DNA barcoding"
    )
    assert meta["journal"] == "Journal of 50% Things"
    assert [a["surname"] for a in meta["authors"]] == ["Müller", "Mańko"]
    for field in ("title", "journal"):
        assert "\\" not in meta[field]


# --- export is the inverse ----------------------------------------------------

@pytest.mark.parametrize("plain,escaped", [
    ("Jammu & Kashmir", r"Jammu \& Kashmir"),
    ("50% of taxa", r"50\% of taxa"),
    ("gene_1", r"gene\_1"),
])
def test_export_escapes_the_specials(plain, escaped):
    """A bare ``&`` breaks the .bib the moment it reaches LaTeX, and a bare
    ``%`` comments out the rest of the line — so an exported file would
    lose fields when any real BibTeX tool re-read it."""
    assert _escape(plain) == escaped


@pytest.mark.parametrize("original", [
    "Lineae chinensis (R.Br.) A.Braun & Vatke",
    "50% of taxa in gene_1 cost $5 (#3)",
    "Müller & Mańko on Nanomia bijuga",
])
def test_round_trip_is_stable(original):
    """#177's acceptance criterion: export → import leaves a served title
    identical to its pre-export value."""
    round_bib = "@article{k, title = {" + _escape(original) + "}, file = {X.pdf}}"
    back = bib_entry_to_metadata(parse_bibtex(round_bib)[0], "X.pdf")["title"]
    assert back == original


def test_literal_braces_do_not_survive_a_round_trip():
    """Known and deliberate: brace removal is total, because bibtex braces
    carry case-protection rather than content, and #141 depends on an
    escaped brace leaving nothing behind. A title with a *literal* brace
    therefore loses it. Documented rather than fixed — the two behaviors
    are in direct conflict and OCR-emitted stray braces are far commoner in
    this corpus than titles that genuinely contain one."""
    original = "A title with {literal} braces"
    round_bib = "@article{k, title = {" + _escape(original) + "}, file = {X.pdf}}"
    back = bib_entry_to_metadata(parse_bibtex(round_bib)[0], "X.pdf")["title"]
    assert back == "A title with literal braces"


# --- #141 must not regress ----------------------------------------------------

def test_escaped_brace_recovery_still_holds():
    """The #141 behavior this function also carries: an escaped brace from
    Grobid OCR leaves no stray backslash."""
    assert unescape_bib_value(r"Des, Ej\{aims") == "Des, Ejaims"
    assert "\\" not in unescape_bib_value(r"a \{b\} c")

"""BibTeX parser robustness — silent truncation and its root cause (#141).

The reported failure: importing a 19,834-entry export parsed only 2,258
entries and reported a clean-looking summary. One unbalanced brace made
the streaming scanner stop, discarding ~88% of the file, and the only
signal was a single WARNING line.

Two distinct defects, tested separately below because they fail
independently:

1. **The imbalance was usually self-inflicted.** Neither brace scanner
   honored backslash escapes, so Grobid OCR output like
   ``author = {Des, Ej\\{aims}`` counted the escaped brace as an opening
   one and the entry never closed. The quoted-string scanner had always
   handled backslashes; only the brace path was missing it.
2. **A genuinely malformed entry aborted the whole parse** instead of
   costing one entry, and the shortfall was not reconciled against the
   number of ``@`` markers in the file.
"""
from __future__ import annotations

import logging

from bib.parser import (
    BibIndex,
    _parse_fields,
    _strip_outer_braces,
    bib_entry_to_metadata,
    count_entry_markers,
    parse_bibtex,
)


# --- 1. escaped braces are not brace-counting events --------------------------


def test_escaped_open_brace_does_not_swallow_the_file():
    """The #141 repro. Before the fix this returned only ``good1``."""
    text = r"""@article{good1,
      author = {Smith, A},
      title = {Fine},
    }
    @article{bad,
      author = {Des, Ej\{aims},
      title = {Broken},
    }
    @article{good2,
      author = {Jones, B},
      title = {Also fine},
    }
    """
    entries = parse_bibtex(text)
    assert [e["_key"] for e in entries] == ["good1", "bad", "good2"]
    # ...and the offending entry parses fully, not just up to the escape.
    bad = next(e for e in entries if e["_key"] == "bad")
    assert bad["title"] == "Broken"


def test_escaped_close_brace_is_literal():
    text = r"""@article{k, author = {A lone \} here}, year = {1999}, }"""
    entries = parse_bibtex(text)
    assert len(entries) == 1
    assert entries[0]["year"] == "1999"  # the field *after* the escape survives


def test_escaped_braces_leave_no_stray_backslash_in_metadata():
    """`_strip_outer_braces` removes braces; an escaped brace must not
    leave its backslash behind in the JSON metadata."""
    assert _strip_outer_braces(r"Des, Ej\{aims") == "Des, Ejaims"
    assert "\\" not in _strip_outer_braces(r"a \{b\} c")


def test_nested_braces_still_work():
    """Protected casing is the reason braces nest at all — don't regress it."""
    text = """@article{k, title = {The {DNA} of {Nanomia bijuga}}, }"""
    entries = parse_bibtex(text)
    assert entries[0]["title"] == "The {DNA} of {Nanomia bijuga}"
    assert _strip_outer_braces(entries[0]["title"]) == "The DNA of Nanomia bijuga"


def test_quoted_values_unaffected():
    text = '''@article{k, title = "A quoted \\" title", year = {2001}, }'''
    entries = parse_bibtex(text)
    assert entries[0]["year"] == "2001"


# --- 2. recovery + reconciliation ---------------------------------------------


def test_genuinely_malformed_entry_costs_one_entry_not_the_file():
    text = """@article{a1, title = {One}, }
@article{a2, title = {Two {unclosed }, }
@article{a3, title = {Three}, }
"""
    entries = parse_bibtex(text)
    keys = [e["_key"] for e in entries]
    # a2 is unrecoverable, but a3 — after it — must still be parsed.
    assert "a1" in keys and "a3" in keys
    assert "a2" not in keys


def test_malformed_entry_is_counted_and_reported(caplog):
    text = """@article{a1, title = {One}, }
@article{a2, title = {Two {unclosed }, }
@article{a3, title = {Three}, }
"""
    with caplog.at_level(logging.WARNING, logger="bib.parser"):
        parse_bibtex(text)
    joined = "\n".join(r.getMessage() for r in caplog.records)
    # The reconciliation line is what makes a shortfall visible: the parsed
    # count alone always looks plausible, which is why #141 went unnoticed.
    assert "skipped as malformed" in joined, joined
    assert "'@' markers in file" in joined, joined


def test_no_warning_on_a_clean_file(caplog):
    """A well-formed bib must stay silent, or the new warnings become noise
    that operators learn to ignore."""
    text = """@article{a1, title = {One}, }
@article{a2, title = {Two}, }
"""
    with caplog.at_level(logging.WARNING, logger="bib.parser"):
        entries = parse_bibtex(text)
    assert len(entries) == 2
    assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []


def test_unbalanced_field_warns_and_names_the_entry(caplog):
    """The quiet sibling: an unbalanced brace inside a *value* swallows the
    entry's remaining fields, and used to do so with no message at all.

    Exercised against `_parse_fields` directly rather than through
    `parse_bibtex`, because the same brace arithmetic governs both levels —
    a body whose field is unbalanced makes the *entry* unbalanced too, so
    the entry-level guard fires first and this path is unreachable from a
    whole-file parse. It is still worth covering: `_parse_fields` is called
    with a body the entry scanner has already accepted, and a future change
    to either scanner could separate the two.
    """
    with caplog.at_level(logging.WARNING, logger="bib.parser"):
        fields = _parse_fields("author = {A {unclosed, title = {Lost}",
                               cite_key="halfparsed")
    joined = "\n".join(r.getMessage() for r in caplog.records)
    assert joined, "an unbalanced field must not be silent"
    assert "halfparsed" in joined, joined
    assert "author" in joined, joined
    # The later field really is lost — that is the damage being reported.
    assert "title" not in fields


def test_count_entry_markers_ignores_non_records():
    text = """@comment{ignored}
@string{x = {y}}
@article{real1, title = {t}, }
@article{real2, title = {t}, }
"""
    assert count_entry_markers(text) == 4      # crude by design: counts all
    assert len(parse_bibtex(text)) == 2       # ...but only records parse


def test_comment_lines_still_skipped():
    text = """% @article{commented_out, title = {no}, }
@article{real, title = {yes}, }
"""
    entries = parse_bibtex(text)
    assert [e["_key"] for e in entries] == ["real"]


# --- end-to-end through the consumer -----------------------------------------


def test_index_and_metadata_survive_a_damaged_file():
    """The whole point: a bib with one bad entry still indexes every good
    entry, so an import no longer silently loses most of a library."""
    text = r"""@article{p1, file = {One.pdf}, title = {First}, year = {2001}, }
@article{ocr_damaged, author = {X, Ej\{aims}, file = {Two.pdf}, title = {Second}, year = {2002}, }
@article{p3, file = {Three.pdf}, title = {Third}, year = {2003}, }
"""
    idx = BibIndex(parse_bibtex(text))
    assert len(idx) == 3
    for name in ("One.pdf", "Two.pdf", "Three.pdf"):
        assert name in idx, f"{name} missing from the index"
    md = bib_entry_to_metadata(idx.lookup("Two.pdf"), "Two.pdf")
    assert md["year"] == 2002
    assert md["title"] == "Second"
    assert "\\" not in md["authors"][0]["surname"]

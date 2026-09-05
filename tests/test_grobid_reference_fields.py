"""TEI bibliography fields keep article/book titles distinct from journals."""
from __future__ import annotations

from pipeline.grobid_client import parse_tei_references


def _tei(entries: str) -> str:
    return f"""\
<TEI xmlns="http://www.tei-c.org/ns/1.0">
  <text><back><listBibl>{entries}</listBibl></back></text>
</TEI>
"""


def test_journal_only_reference_does_not_copy_journal_into_work_title():
    refs = parse_tei_references(_tei("""
      <biblStruct xml:id="b0">
        <monogr>
          <title level="j">Phytochemistry</title>
          <imprint><date when="1989"/></imprint>
        </monogr>
      </biblStruct>
    """))

    assert refs[0]["title"] == ""
    assert refs[0]["journal"] == "Phytochemistry"


def test_article_title_wins_and_journal_remains_separate():
    refs = parse_tei_references(_tei("""
      <biblStruct xml:id="b0">
        <analytic><title level="a">Iridoid glucosides from Viburnum</title></analytic>
        <monogr><title level="j">Phytochemistry</title></monogr>
      </biblStruct>
    """))

    assert refs[0]["title"] == "Iridoid glucosides from Viburnum"
    assert refs[0]["journal"] == "Phytochemistry"


def test_non_journal_monograph_title_is_a_valid_book_title():
    refs = parse_tei_references(_tei("""
      <biblStruct xml:id="b0">
        <monogr><title level="m">Flora of Turkey</title></monogr>
      </biblStruct>
    """))

    assert refs[0]["title"] == "Flora of Turkey"
    assert refs[0]["journal"] == ""


def test_unlevelled_duplicate_journal_is_not_a_title():
    refs = parse_tei_references(_tei("""
      <biblStruct xml:id="b0">
        <monogr>
          <title>Nature</title>
          <title level="j">Nature</title>
        </monogr>
      </biblStruct>
    """))

    assert refs[0]["title"] == ""
    assert refs[0]["journal"] == "Nature"


def test_raw_citation_is_preserved_independently_of_lossy_structured_fields():
    refs = parse_tei_references(_tei("""
      <biblStruct xml:id="b50">
        <analytic><title level="a">Diversity and vertical distribution</title>
          <author><persName><forename>K</forename><surname>Hopcroft</surname></persName></author>
        </analytic>
        <note type="raw_reference">Kosobokova, K., Hopcroft, R.R., 2010.
          Diversity and vertical distribution.</note>
      </biblStruct>
    """))
    assert refs[0]["authors"] == ["K Hopcroft"]
    assert refs[0]["raw"] == ("Kosobokova, K., Hopcroft, R.R., 2010. "
                              "Diversity and vertical distribution.")

"""Expansion of abbreviated genus binomials — `Ph. pelagica` (#164).

Taxonomic literature abbreviates the genus after first mention, so for a
corpus of original descriptions this is the central case rather than an
edge one: the paper that *erects* a species is the one least likely to
spell the genus out on every line. Olfers 1824 is a five-species key for
*Physalia* that yielded one genus-level taxon and no species at all.

Naive expansion would be worse than the under-extraction it replaces,
because `Ph.` is genuinely ambiguous here — *Physalia* and *Physophora*
are both in this corpus. Three gates keep it honest: the genus must be
written out in full somewhere in the same document, the expansion must
be a name in the taxonomy snapshot, and the survivors must agree on one
accepted taxon.
"""
from __future__ import annotations

import sqlite3

import pytest

from pipeline.taxa import TaxonomyDB, extract_taxon_mentions
from pipeline.taxonomy_ingest import create_schema, insert_records, make_record


@pytest.fixture
def taxonomy(tmp_path):
    """A snapshot with the ambiguity the issue warns about, and the
    historical-spelling shape that made the first cut over-report."""
    path = tmp_path / "taxonomy.sqlite"
    conn = sqlite3.connect(path)
    create_schema(conn)
    insert_records(conn, [
        # `Physalis` is a 19th-century spelling of the genus and lives in
        # the snapshot as a name of its own, exactly as in WoRMS. Olfers
        # 1824 prints both.
        make_record(taxon_id="1", scientific_name="Physalia",
                    taxon_rank="Genus", taxonomic_status="accepted",
                    extra_names=["Physalis"]),
        make_record(taxon_id="2", scientific_name="Physalia physalis",
                    taxonomic_status="accepted",
                    # Two historical spellings of one taxon, which is not
                    # ambiguity — plus a real synonym epithet.
                    extra_names=["Physalia pelagica", "Physalis pelagica"]),
        make_record(taxon_id="3", scientific_name="Physophora",
                    taxon_rank="Genus", taxonomic_status="accepted"),
        make_record(taxon_id="4", scientific_name="Physophora hydrostatica",
                    taxonomic_status="accepted"),
        # `Ph. megalista` is the genuinely ambiguous case: an accepted
        # species under one spelling, a synonym of another taxon under
        # the other.
        make_record(taxon_id="5", scientific_name="Physalia megalista",
                    taxonomic_status="accepted"),
        make_record(taxon_id="6", scientific_name="Physalis megalista",
                    taxonomic_status="accepted"),
    ])
    conn.commit()
    conn.close()
    return TaxonomyDB(path)


def _chunks(*texts):
    return [{"chunk_id": f"chunk_{i}", "text": t} for i, t in enumerate(texts)]


def test_an_abbreviation_expands_against_a_genus_named_in_full(taxonomy):
    out = extract_taxon_mentions(
        _chunks("Physalia Arethusa Til. 2. Ph. pelagica. Subovalis, altera"),
        taxonomy,
    )
    expanded = [m for m in out["mentions"] if m.get("method") == "abbreviated_genus"]
    assert len(expanded) == 1
    assert out["abbreviations_expanded"] == 1
    assert expanded[0]["accepted_name"] == "Physalia physalis"
    # The printed form and the resolved name are kept apart, so nothing
    # downstream takes an inferred name for an observed one.
    assert expanded[0]["mention_text"] == "Ph. pelagica"
    assert expanded[0]["matched_text"] == "Physalia pelagica"
    assert expanded[0]["expanded_from"] == "Ph"


def test_an_ocr_comma_for_the_period_still_expands(taxonomy):
    """Olfers 1824 prints `Ph, pelagica` — the period lost to OCR."""
    out = extract_taxon_mentions(
        _chunks("Physalia Arethusa. 2. Ph, pelagica. Subovalis"), taxonomy)
    assert out["abbreviations_expanded"] == 1


def test_nothing_is_expanded_without_the_genus_in_the_document(taxonomy):
    """The whole safeguard: an abbreviation resolves against this
    document's own genera, not against the taxonomy at large."""
    out = extract_taxon_mentions(_chunks("2. Ph. pelagica. Subovalis"), taxonomy)
    assert out["abbreviations_expanded"] == 0
    assert out["total_mentions"] == 0


def test_a_genuine_ambiguity_records_nothing_but_reports_itself(taxonomy):
    """Two spellings that resolve to *different* accepted taxa. Writing a
    confident guess into taxon_mentions.sqlite is worse than the
    under-extraction it replaces, so nothing is written — but the
    ambiguity is reported rather than dropped in silence."""
    out = extract_taxon_mentions(
        _chunks("Physalia and Physalis both. 3. Ph. megalista Peron et Lesueur."),
        taxonomy,
    )
    assert out["abbreviations_expanded"] == 0
    assert out["abbreviations_unresolved_count"] == 1
    report = out["abbreviations_unresolved"][0]
    assert report["mention_text"] == "Ph. megalista"
    assert report["reason"] == "ambiguous_abbreviation"
    assert sorted(report["candidates"]) == ["Physalia megalista",
                                            "Physalis megalista"]


def test_two_spellings_of_one_taxon_are_not_an_ambiguity(taxonomy):
    """Historical spellings live in the snapshot as names of their own.
    Olfers 1824 prints both *Physalia* and *Physalis*, and `Physalia
    pelagica` and `Physalis pelagica` are two names for one taxon —
    counting name strings called that ambiguous and dropped a mention
    that was never in doubt."""
    out = extract_taxon_mentions(
        _chunks("Physalia here. Physalia again. Physalis once. 2. Ph. pelagica."),
        taxonomy,
    )
    assert out["abbreviations_expanded"] == 1
    assert out["abbreviations_unresolved_count"] == 0
    expanded = [m for m in out["mentions"] if m.get("method") == "abbreviated_genus"][0]
    # The document's own preference picks the representative spelling.
    assert expanded["matched_text"] == "Physalia pelagica"


def test_an_ambiguous_prefix_is_split_by_the_epithet(taxonomy):
    """`Ph.` matches both genera in this corpus, so document context
    alone would be guessing — but the taxonomy knows `Physophora
    hydrostatica` is a name and `Physalia hydrostatica` is not."""
    out = extract_taxon_mentions(
        _chunks("Physalia and Physophora are both here. Ph. hydrostatica is it."),
        taxonomy,
    )
    expanded = [m for m in out["mentions"] if m.get("method") == "abbreviated_genus"]
    assert len(expanded) == 1
    assert expanded[0]["matched_text"] == "Physophora hydrostatica"


def test_a_bibliographic_abbreviation_is_not_a_taxon(taxonomy):
    """The regex is loose on purpose — `Taf. iii` matches it. The gates
    behind it are what decide, and a plate reference clears none."""
    out = extract_taxon_mentions(
        _chunks("Physalia physalis, Taf. iii and No. seven, cf. supra"), taxonomy)
    assert out["abbreviations_expanded"] == 0


def test_an_expansion_does_not_overlap_a_name_written_out_in_full(taxonomy):
    """`Physalia physalis` must be read once, as itself, not also as an
    abbreviation hiding inside it."""
    out = extract_taxon_mentions(_chunks("Physalia physalis Linnaeus"), taxonomy)
    assert out["abbreviations_expanded"] == 0
    assert out["total_mentions"] == 1


def test_expansions_carry_their_own_method_for_downstream_filtering(taxonomy):
    """An expansion is an inference from the document's genus list, not a
    name read off the page. A query that needs verbatim evidence has to be
    able to tell them apart."""
    out = extract_taxon_mentions(
        _chunks("Physalia Arethusa. 2. Ph. pelagica."), taxonomy)
    methods = {m.get("method", "regex_taxonomy") for m in out["mentions"]}
    assert methods == {"regex_taxonomy", "abbreviated_genus"}


def test_mentions_stay_in_text_order_within_a_chunk(taxonomy):
    """Pass 2 interleaves rather than appending a block at the end."""
    out = extract_taxon_mentions(
        _chunks("Ph. pelagica first, then Physalia spelled out, "
                "then Ph. pelagica again"),
        taxonomy,
    )
    spans = [m["text_span"][0] for m in out["mentions"]]
    assert spans == sorted(spans)
    assert out["abbreviations_expanded"] == 2


def test_the_unresolved_report_is_bounded(taxonomy):
    """taxa.json ships in the served bundle, so the ambiguity report is a
    diagnostic, not a parallel mention list."""
    text = "Physalia and Physalis. " + "Ph. megalista. " * 80
    out = extract_taxon_mentions(_chunks(text), taxonomy)
    assert out["abbreviations_unresolved_count"] == 80
    assert len(out["abbreviations_unresolved"]) == 50

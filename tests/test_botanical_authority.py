"""Botanical (ICN) authorship, and reporting a capability that cannot work.

Authority linking matches a taxon's authorship against a work by author
*and year*, which is the zoological (ICZN) convention WoRMS supplies.
Botanical authorship is author-only by correct citation practice —
`Rehder`, `(Kache) Hesse` — so there is no year to match on. On a
702-paper *Viburnum* corpuscle that meant 0 links from 889 authorship
strings, and `get_original_description` answering `null` for every taxon
while presenting itself as an available tool (#175).

The author-only matching path the issue proposes was measured on that
corpuscle before being declined: pairing the authority surname with the
epithet appearing in a work title yields 3 candidates from 889 taxa, and
2 of the 3 are wrong. So this fixes the reporting, deliberately, and the
decision is recorded in `_record_authority_convention`.
"""
import json
import sqlite3

import pytest

from bib import authority as a


# --- parsing ------------------------------------------------------------


@pytest.mark.parametrize("authority,surnames", [
    # Real values from the World Checklist of Vascular Plants snapshot.
    ("Rehder", ["Rehder"]),
    ("Rouleau", ["Rouleau"]),
    ("House", ["House"]),
    ("Nakai", ["Nakai"]),
    # The parenthesised author published the original name; the one after
    # it published this combination. Both are returned, in printed order.
    ("(Rehder) Rehder", ["Rehder"]),
    ("(Kache) Hesse", ["Kache", "Hesse"]),
    # Botanical initials run together without spaces, which the ICZN
    # path's initial-stripping never had to handle.
    ("(Huxley) P.S.Hsu", ["Huxley", "Hsu"]),
    ("(Hiyama) H.Hara", ["Hiyama", "Hara"]),
    ("(A.Juss.) Nakai", ["Juss", "Nakai"]),
    # `Vent.` is an abbreviated surname (Ventenat), not an initial.
    ("(Vent.) P.Silva", ["Vent", "Silva"]),
    # `L.` for Linnaeus is a whole authorship.
    ("L.", ["L"]),
])
def test_icn_author_only_authorship_parses(authority, surnames):
    parsed = a.parse_authority(authority)
    assert parsed is not None, authority
    assert parsed == (surnames, None)


@pytest.mark.parametrize("authority,surnames,year", [
    ("Eschscholtz, 1829", ["Eschscholtz"], 1829),
    ("(Huxley, 1859)", ["Huxley"], 1859),
    ("Quoy & Gaimard, 1833", ["Quoy", "Gaimard"], 1833),
    ("L. Agassiz, 1862", ["Agassiz"], 1862),
    ("Lens & van Riemsdijk, 1908", ["Lens", "van Riemsdijk"], 1908),
])
def test_the_zoological_shape_is_unchanged(authority, surnames, year):
    assert a.parse_authority(authority) == (surnames, year)


@pytest.mark.parametrize("authority", ["", None, "???", "--", "1829"])
def test_nothing_name_like_is_still_unparseable(authority):
    """A bare year is a malformed ICZN string, not an ICN one — do not
    silently reinterpret it as an author called "1829"."""
    assert a.parse_authority(authority) is None


# --- the recorded verdict ----------------------------------------------


def _taxonomy(path, rows):
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE IF NOT EXISTS taxa"
                 "(taxon_id TEXT, scientific_name_authorship TEXT)")
    conn.executemany("INSERT INTO taxa VALUES (?,?)", rows)
    conn.commit()
    conn.close()


def _db():
    conn = sqlite3.connect(":memory:")
    a.create_schema(conn)
    return conn


def _verdict(conn):
    row = conn.execute(
        "SELECT value FROM build_meta WHERE key=?",
        (a.AUTHORITY_CONVENTION_KEY,),
    ).fetchone()
    return json.loads(row[0]) if row and row[0] else None


def test_a_botanical_taxonomy_is_recorded_as_unsupported(tmp_path, caplog):
    tx = tmp_path / "taxonomy.sqlite"
    _taxonomy(tx, [("1", "Rehder"), ("2", "(Kache) Hesse"), ("3", "Nakai")])
    conn = _db()
    with caplog.at_level("WARNING"):
        a.phase3_authority_links(conn, tx)

    verdict = _verdict(conn)
    assert verdict["supported"] is False
    assert verdict["convention"] == "botanical"
    assert verdict["author_only"] == 3
    assert verdict["year_bearing"] == 0
    assert verdict["links"] == 0
    # A build whose taxonomy can never link has to say so once, loudly.
    assert "Authority linking is unavailable" in caplog.text
    assert "#175" in caplog.text


def test_a_zoological_taxonomy_is_recorded_as_supported(tmp_path, caplog):
    tx = tmp_path / "taxonomy.sqlite"
    _taxonomy(tx, [("1", "Eschscholtz, 1829"), ("2", "(Huxley, 1859)")])
    conn = _db()
    with caplog.at_level("WARNING"):
        a.phase3_authority_links(conn, tx)

    verdict = _verdict(conn)
    assert verdict["supported"] is True
    assert verdict["convention"] == "zoological"
    assert verdict["year_bearing"] == 2
    assert "Authority linking is unavailable" not in caplog.text


def test_a_yearless_authorship_creates_no_stub_work(tmp_path):
    """The point of declining the author-only path: a botanical taxonomy
    must not fill `works` with yearless stubs, or `taxon_work_links` with
    guesses."""
    tx = tmp_path / "taxonomy.sqlite"
    _taxonomy(tx, [("1", "Rehder"), ("2", "(Vent.) P.Silva")])
    conn = _db()
    a.phase3_authority_links(conn, tx)
    assert conn.execute("SELECT COUNT(*) FROM works").fetchone()[0] == 0
    assert conn.execute(
        "SELECT COUNT(*) FROM taxon_work_links").fetchone()[0] == 0


def test_a_mixed_taxonomy_links_what_it_can_and_reports_supported(tmp_path):
    """One year-bearing string is enough for the capability to exist."""
    tx = tmp_path / "taxonomy.sqlite"
    _taxonomy(tx, [("1", "Rehder"), ("2", "Eschscholtz, 1829")])
    conn = _db()
    a.phase3_authority_links(conn, tx)
    verdict = _verdict(conn)
    assert verdict["supported"] is True
    assert (verdict["year_bearing"], verdict["author_only"]) == (1, 1)


def test_recording_the_verdict_does_not_break_the_no_op_guarantee(tmp_path):
    """Phase 3 promises a no-op run performs no writes, and the verdict
    is written on every run — so it has to be read before it is written."""
    tx = tmp_path / "taxonomy.sqlite"
    _taxonomy(tx, [("1", "Rehder")])
    conn = _db()
    a.phase3_authority_links(conn, tx)
    before = conn.total_changes
    assert a.phase3_authority_links(conn, tx) == 0
    assert conn.total_changes == before


# --- what the served tool reports --------------------------------------


def test_the_bundle_accessor_reads_the_verdict_back(tmp_path):
    from mcpsrv.indexes import BiblioAuthority

    db = tmp_path / "biblio_authority.sqlite"
    tx = tmp_path / "taxonomy.sqlite"
    _taxonomy(tx, [("1", "Rehder")])
    conn = sqlite3.connect(db)
    a.create_schema(conn)
    a.phase3_authority_links(conn, tx)
    conn.close()

    assert BiblioAuthority(db).authority_linking()["supported"] is False


def test_a_bundle_predating_the_verdict_reports_unknown(tmp_path):
    """Absence of the row is not evidence either way, and must not be
    reported as "unsupported" on an older bundle."""
    from mcpsrv.indexes import BiblioAuthority

    db = tmp_path / "biblio_authority.sqlite"
    conn = sqlite3.connect(db)
    a.create_schema(conn)
    conn.commit()
    conn.close()

    linking = BiblioAuthority(db).authority_linking()
    assert linking["supported"] is None
    assert linking["convention"] == "unknown"


# --- what the served tool reports, end to end --------------------------


@pytest.fixture
def serving(request):
    """Install a minimal index so get_original_description is callable."""
    import types

    from mcpsrv import app as mcp_app

    linking, works = request.param

    def _install():
        taxonomy_db = types.SimpleNamespace(
            lookup=lambda name: {
                "accepted_taxon_id": "t:1",
                "accepted_name": "Viburnum dentatum",
                "rank": "species",
            } if name else None
        )
        biblio_db = types.SimpleNamespace(
            work_for_taxon=lambda _tid: [dict(w) for w in works],
            authority_linking=lambda: dict(linking),
            get_authors=lambda _wid: [],
            citation_count=lambda _wid: 0,
        )
        return types.SimpleNamespace(taxonomy_db=taxonomy_db,
                                     biblio_db=biblio_db)

    original = mcp_app._INDEX
    mcp_app.set_index(_install())
    yield
    mcp_app.set_index(original)


_BOTANICAL = ({"supported": False, "convention": "botanical",
               "author_only": 889, "authorship_strings": 889}, [])
_ZOOLOGICAL_MISS = ({"supported": True, "convention": "zoological",
                     "author_only": 0, "authorship_strings": 800}, [])
_UNKNOWN = ({"supported": None, "convention": "unknown"}, [])


@pytest.mark.parametrize("serving", [_BOTANICAL], indirect=True)
def test_the_tool_reports_an_unsupported_capability(serving):
    """The whole point of #175: `null` with "no matching work found" reads
    as "no such paper exists", which is a different claim entirely."""
    from mcpsrv.tools.bibliography import get_original_description

    out = get_original_description("Viburnum dentatum")
    assert out["original_description"] is None
    assert out["unsupported"] is True
    assert out["reason_code"] == "authority_convention_unsupported"
    assert "botanical" in out["note"]
    assert "not evidence" in out["note"]
    assert out["authority_linking"]["author_only"] == 889


@pytest.mark.parametrize("serving", [_ZOOLOGICAL_MISS], indirect=True)
def test_a_genuine_miss_still_reads_as_a_miss(serving):
    """Where linking does work, an empty answer means what it says."""
    from mcpsrv.tools.bibliography import get_original_description

    out = get_original_description("Marrus claudanielis")
    assert out["original_description"] is None
    assert "unsupported" not in out
    assert out["note"] == "no matching work found in the authority database"


@pytest.mark.parametrize("serving", [_UNKNOWN], indirect=True)
def test_an_older_bundle_is_not_declared_unsupported(serving):
    """A bundle built before the verdict existed has no row; that is not
    evidence the capability is missing."""
    from mcpsrv.tools.bibliography import get_original_description

    out = get_original_description("Marrus claudanielis")
    assert "unsupported" not in out
    assert out["authority_linking"]["convention"] == "unknown"

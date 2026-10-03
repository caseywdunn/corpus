"""Build-time evidence for original-description candidates (#311).

Author/year agreement supplies a lead, never proof of an original description.
Opening-text markers can support a lead without converting it into a curator
verdict. The server reads these bounded records; it does not scan documents.
"""
from __future__ import annotations

import hashlib
import json
import re

PRODUCER = "taxon-authority-candidates-v2"


def create_schema(conn):
    conn.execute("""CREATE TABLE IF NOT EXISTS taxon_authority_candidates (
        taxon_id TEXT NOT NULL,
        work_id TEXT NOT NULL REFERENCES works(work_id),
        confidence REAL NOT NULL,
        basis_json TEXT NOT NULL,
        producer_version TEXT NOT NULL,
        PRIMARY KEY(taxon_id,work_id)
    )""")


def description_evidence(output_dir, corpus_hash, scientific_name):
    if not corpus_hash or not scientific_name:
        return None
    names = scientific_name.split()
    if len(names) < 2:
        return None
    # Full binomial or a printed genus initial, followed immediately by an
    # explicit new-species marker. A taxon mention alone is not evidence.
    genus, epithet = names[:2]
    taxon = rf"(?:{re.escape(genus)}|{re.escape(genus[0])}\s*\.)\s+{re.escape(epithet)}"
    pattern = re.compile(rf"\b{taxon}\s+(?:sp\s*\.\s*nov\s*\.?|new\s+species)\b", re.IGNORECASE)
    path = output_dir / "documents" / corpus_hash / "text.json"
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    text = payload.get("text")
    if not isinstance(text, str):
        return None
    match = pattern.search(text[:6000])
    if not match:
        return None
    return {"kind": "opening_text_new_species_marker", "corpus_hash": corpus_hash,
            "source": "text.json", "text_sha256": hashlib.sha256(text.encode()).hexdigest(),
            "char_start": match.start(), "char_end": match.end(),
            "excerpt": text[max(0, match.start()-90):min(len(text), match.end()+90)]}


def candidates_for(conn, surnames, year, *, output_dir, scientific_name):
    from .authority import normalize_for_key
    expected = [normalize_for_key(s) for s in surnames]
    rows = conn.execute("""SELECT DISTINCT w.work_id,w.corpus_hash,w.title
        FROM works w JOIN work_authors a ON a.work_id=w.work_id
        WHERE a.position=0 AND a.surname_normalized=? AND w.year=?
          AND COALESCE(TRIM(w.title),'') != ''
          AND (w.source!='taxon_authority' OR w.in_corpus=1 OR w.bib_imported_at IS NOT NULL)
        ORDER BY w.work_id""", (expected[0], year)).fetchall()
    result = {}
    for work_id, corpus_hash, _title in rows:
        actual = [r[0] for r in conn.execute("SELECT surname_normalized FROM work_authors WHERE work_id=? ORDER BY position", (work_id,))]
        exact = actual == expected
        # First-author/year candidates remain visible when later author
        # parsing differs, but the disagreement is explicit and lowers rank.
        basis = {"kind": "ordered_authors_year" if exact else "first_author_year",
                 "authority_authors": surnames, "year": year,
                 "complete_author_list_match": exact,
                 "requires_source_review": True}
        confidence = 0.85 if exact else 0.6
        evidence = description_evidence(output_dir, corpus_hash, scientific_name)
        if evidence:
            basis["source_evidence"] = evidence
            confidence = 0.95 if exact else 0.75
        result[work_id] = (confidence, json.dumps(basis, sort_keys=True, ensure_ascii=False), PRODUCER)
    return result


def replace_current(conn, desired):
    current = {(taxon_id, work_id): (confidence, basis, producer) for taxon_id, work_id, confidence, basis, producer in conn.execute(
        "SELECT taxon_id,work_id,confidence,basis_json,producer_version FROM taxon_authority_candidates")}
    for key in current.keys()-desired.keys():
        conn.execute("DELETE FROM taxon_authority_candidates WHERE taxon_id=? AND work_id=?", key)
    for key, value in desired.items():
        if current.get(key) != value:
            conn.execute("INSERT OR REPLACE INTO taxon_authority_candidates VALUES (?,?,?,?,?)", (*key, *value))
    return len(current.keys()-desired.keys()) + sum(current.get(k) != v for k,v in desired.items())

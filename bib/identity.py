"""Citation-unit identity, distinct from shared DOI membership (#299, #300)."""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict

DOCUMENT_IDENTITY_PRODUCER = "document-identity-v2"


def normalized_identity(meta):
    from .authority import normalize_for_key
    fields = ("title", "year", "journal", "volume", "number", "pages", "eid", "articleno", "booktitle", "chapter", "edition")
    result = {key: normalize_for_key(str(meta.get(key) or "")) for key in fields}
    # A DOI and complete title can establish one work despite an extracted
    # author's typo. Without a DOI the full author list helps disambiguate
    # same-title series parts; forename spelling is field authority, not ID.
    if not meta.get("doi"):
        result["authors"] = [normalize_for_key(a.get("surname", "")) for a in meta.get("authors", [])]
    return result


def compatible(left, right):
    a, b = normalized_identity(left), normalized_identity(right)
    return all(not a.get(key) or not b.get(key) or a[key] == b[key] for key in a.keys() | b.keys())


def base_identity(meta, corpus_hash):
    from .authority import make_corpus_guid, normalize_doi
    doi = normalize_doi(meta.get("doi") or "")
    surname = next((a.get("surname") for a in meta.get("authors", []) if a.get("surname")), "")
    if doi:
        return doi, "doi"
    if surname:
        return make_corpus_guid(surname, meta.get("year"), meta.get("title") or meta.get("filename") or corpus_hash), "corpus_key"
    return f"corpus:{corpus_hash}", "corpus_key"


def document_identities(inputs):
    """Plan all IDs before ingest so collisions cannot depend on traversal order.

    Unknown fields do not split otherwise compatible scans. A partially
    described record compatible with multiple contradictory groups gets its
    own identity: missing evidence cannot choose a part for the caller.
    """
    groups = defaultdict(list)
    for folder, meta in inputs:
        base, kind = base_identity(meta, folder.name)
        groups[(base, kind)].append((folder.name, meta))
    result = {}
    for (base, kind), members in groups.items():
        clusters = []
        # More complete evidence first, then content; never ingestion order.
        ordered = sorted(members, key=lambda pair: (
            -sum(bool(v) for v in normalized_identity(pair[1]).values()),
            json.dumps(normalized_identity(pair[1]), sort_keys=True), pair[0]))
        for sha, meta in ordered:
            matches = [cluster for cluster in clusters if all(compatible(meta, item[1]) for item in cluster)]
            if len(matches) == 1:
                matches[0].append((sha, meta))
            else:
                clusters.append([(sha, meta)])
        for cluster in clusters:
            if len(clusters) == 1:
                work_id = base
            else:
                canonical = normalized_identity(cluster[0][1])
                token = hashlib.sha256(json.dumps(canonical, sort_keys=True, ensure_ascii=False).encode()).hexdigest()[:16]
                work_id = base + "#part:" + token
            for sha, _ in cluster:
                result[sha] = (work_id, kind, base, len(clusters) > 1)
    return result


def work_metadata(conn, work_id):
    cur = conn.execute("SELECT * FROM works WHERE work_id=?", (work_id,))
    row = cur.fetchone()
    if not row:
        return {}
    result = dict(zip([d[0] for d in cur.description], row))
    result["authors"] = [{"surname": s, "forename": f or ""} for s, f in conn.execute(
        "SELECT surname,forename FROM work_authors WHERE work_id=? ORDER BY position", (work_id,))]
    return result


def select_reference_candidate(conn, work_ids, ref):
    """Choose one part from shared-key candidates using positive metadata evidence."""
    from .authority import extract_surname_from_ref_author
    meta = dict(ref, authors=[{"surname": extract_surname_from_ref_author(a)} for a in ref.get("authors", [])])
    rows = [(wid, work_metadata(conn, wid)) for wid in work_ids]
    candidates = [(wid, candidate) for wid, candidate in rows if compatible(candidate, meta)]
    local = [wid for wid, candidate in candidates if candidate.get("in_corpus")]
    if local:
        return local[0] if len(local) == 1 else None
    return candidates[0][0] if len(candidates) == 1 else None


def record_decision(conn, corpus_hash, source_id, target_id, reason, evidence):
    conn.execute("""CREATE TABLE IF NOT EXISTS work_identity_decisions (
        corpus_hash TEXT NOT NULL, source_work_id TEXT NOT NULL,
        target_work_id TEXT NOT NULL, reason TEXT NOT NULL,
        producer_version TEXT NOT NULL, evidence_json TEXT NOT NULL,
        PRIMARY KEY(corpus_hash,source_work_id,target_work_id,reason,producer_version)
    )""")
    conn.execute("INSERT OR REPLACE INTO work_identity_decisions VALUES (?,?,?,?,?,?)",
                 (corpus_hash, source_id, target_id, reason, DOCUMENT_IDENTITY_PRODUCER,
                  json.dumps(evidence, sort_keys=True, ensure_ascii=False)))

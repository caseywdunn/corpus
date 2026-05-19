"""Bibliographic MCP tools.

Surfaces queries over the corpus-wide ``biblio_authority.sqlite``:
get_bibliography, get_intext_citations, get_excerpts_citing,
get_citation_graph, resolve_reference, get_missing_references,
get_original_description, get_works_by_author.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from ..app import _load_json, _need_index, mcp

if TYPE_CHECKING:
    from ..indexes import BiblioAuthority  # noqa: F401  — annotation only


@mcp.tool()
def get_bibliography(
    paper_hash: str,
    resolved: bool = False,
    limit: int = 50,
    offset: int = 0,
) -> List[Dict]:
    """Parsed references for one paper (from Grobid TEI). Defaults to
    50 references; paginate via ``offset``.

    With ``resolved=True``, each reference is enriched with its
    bibliographic authority record: ``work_id``, ``in_corpus``,
    ``corpus_hash`` (if in corpus), and ``cited_by_count``. Click
    through to a cited work via ``get_paper(corpus_hash)``. Requires
    the bibliographic authority database; falls back to unresolved
    output if unavailable.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [{"error": f"no such paper_hash: {paper_hash}"}]
    refs = _load_json(Path(p["hash_dir"]) / "references.json", default={}) or {}
    ref_list = (refs.get("references", []) or [])[offset: offset + int(limit)]
    if not resolved or idx.biblio_db is None:
        return ref_list

    # Enrich each reference with authority DB info
    for ref in ref_list:
        cur = idx.biblio_db.conn.execute(
            """SELECT c.cited_work_id, w.in_corpus, w.corpus_hash,
                      w.title AS resolved_title, w.doi AS resolved_doi,
                      c.match_method, c.match_score
               FROM citations c JOIN works w ON c.cited_work_id = w.work_id
               WHERE c.citing_corpus_hash = ? AND c.grobid_xml_id = ?""",
            (paper_hash, ref.get("xml_id", "")),
        )
        row = cur.fetchone()
        if row:
            ref["work_id"] = row["cited_work_id"]
            ref["in_corpus"] = bool(row["in_corpus"])
            ref["corpus_hash"] = row["corpus_hash"]
            ref["match_method"] = row["match_method"]
            ref["cited_by_count"] = idx.biblio_db.citation_count(row["cited_work_id"])
        else:
            ref["work_id"] = None
            ref["in_corpus"] = False
            ref["corpus_hash"] = None
    return ref_list


@mcp.tool()
def get_intext_citations(
    paper_hash: str,
    limit: int = 100,
    offset: int = 0,
) -> Dict:
    """In-text citations parsed from this paper's body.

    Returns ``{paragraphs, citations, total_paragraphs,
    total_citations}``. Each citation record carries
    ``target_xml_id`` (links to ``get_bibliography(xml_id=...)``;
    ``None`` for unresolved), ``surface`` (the cited text), the
    enclosing ``section`` heading, and ``para_index`` into the
    returned ``paragraphs`` list.

    ``limit`` + ``offset`` paginate over the citation list. Returned
    ``paragraphs`` is a sublist containing only the paragraphs referenced
    by the citations on this page (deduplicated, in source order); each
    citation's ``para_index`` is remapped to its position in that
    sublist. Paginating over paragraphs and citations independently
    breaks the index relationship, so pagination is anchored to the
    citation list and the paragraph sublist follows.

    Returns ``{"error": ...}`` if intext_citations.json is missing
    (backfill with ``backfill_intext_citations.py``).
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    data = _load_json(Path(p["hash_dir"]) / "intext_citations.json", default=None)
    if data is None:
        return {
            "error": "intext_citations.json missing — backfill with "
                     "`python backfill_intext_citations.py <output_dir>`",
        }
    paragraphs_full = data.get("paragraphs") or []
    citations_full = data.get("citations") or []

    start = int(offset)
    end = start + int(limit)
    citations_page = citations_full[start:end]

    referenced_indices = sorted({
        c.get("para_index")
        for c in citations_page
        if isinstance(c.get("para_index"), int)
        and 0 <= c["para_index"] < len(paragraphs_full)
    })
    global_to_local = {gi: li for li, gi in enumerate(referenced_indices)}
    paragraphs_page = [paragraphs_full[gi] for gi in referenced_indices]

    remapped_citations: List[Dict] = []
    for c in citations_page:
        gi = c.get("para_index")
        c_out = dict(c)
        if isinstance(gi, int) and gi in global_to_local:
            c_out["para_index"] = global_to_local[gi]
        else:
            c_out["para_index"] = None
        remapped_citations.append(c_out)

    return {
        "paragraphs": paragraphs_page,
        "citations": remapped_citations,
        "total_paragraphs": len(paragraphs_full),
        "total_citations": len(citations_full),
    }


@mcp.tool()
def get_excerpts_citing(work_id: str, limit: int = 50) -> Dict:
    """Cross-corpus passages citing ``work_id`` — actual paragraph
    text around each citation marker. Answers "show me every passage
    where <Author Year> is cited" (issue #7).

    ``work_id`` is a DOI / ``corpus:...`` / ``bhl:...`` key — same
    space as ``get_citation_graph``. Only resolved citations
    contribute (~60% of in-text refs); Grobid-unmatched refs don't
    surface here.

    Returns ``{work_id, n_excerpts, excerpts: [{citing_paper_hash,
    citing_paper_title, surface, section, paragraph}, ...]}``.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return {"error": "biblio_authority DB not loaded"}

    rows = idx.biblio_db.conn.execute(
        """SELECT c.citing_corpus_hash, c.grobid_xml_id
           FROM citations c
           WHERE c.cited_work_id = ?""",
        (work_id,),
    ).fetchall()

    excerpts: List[Dict] = []
    for row in rows:
        if len(excerpts) >= limit:
            break
        citing_hash = row["citing_corpus_hash"]
        grobid_id = row["grobid_xml_id"]
        if not citing_hash or not grobid_id:
            continue
        p = idx.papers.get(citing_hash)
        if not p:
            continue
        data = _load_json(
            Path(p["hash_dir"]) / "intext_citations.json", default=None,
        )
        if data is None:
            continue
        paras = data.get("paragraphs", [])
        target_match = "#" + grobid_id
        for c in data.get("citations", []):
            if c.get("target_xml_id") != target_match:
                continue
            pi = c.get("para_index")
            paragraph = paras[pi] if isinstance(pi, int) and 0 <= pi < len(paras) else ""
            excerpts.append({
                "citing_paper_hash": citing_hash,
                "citing_paper_title": p.get("title") or "",
                "surface": c.get("surface", ""),
                "section": c.get("section", ""),
                "paragraph": paragraph,
            })
            if len(excerpts) >= limit:
                break

    return {
        "work_id": work_id,
        "n_excerpts": len(excerpts),
        "excerpts": excerpts,
    }


@mcp.tool()
def get_citation_graph(
    work_id: Optional[str] = None,
    paper_hash: Optional[str] = None,
    direction: str = "both",
    depth: int = 1,
) -> Dict:
    """Citation graph around a work.

    Accepts either a ``work_id`` (DOI, ``bhl:…``, or ``corpus:…``) or
    a ``paper_hash`` (for corpus papers). ``direction`` controls which
    edges to follow: ``"citing"`` (papers that cite this work),
    ``"cited_by"`` (papers cited by this work), or ``"both"``.
    ``depth > 1`` follows transitive citations (max 3).

    Returns the root work plus the citation edges.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return {"error": "bibliographic authority database not configured"}
    # Resolve paper_hash to work_id if needed
    if not work_id and paper_hash:
        w = idx.biblio_db.get_work_by_corpus_hash(paper_hash)
        if not w:
            return {"error": f"no work found for paper_hash: {paper_hash}"}
        work_id = w["work_id"]
    if not work_id:
        return {"error": "provide either work_id or paper_hash"}

    root = idx.biblio_db.get_work(work_id)
    if not root:
        return {"error": f"no such work_id: {work_id}"}

    depth = min(int(depth), 3)
    result: Dict[str, Any] = {
        "root": {
            **root,
            "authors": idx.biblio_db.get_authors(work_id),
            "cited_by_count": idx.biblio_db.citation_count(work_id),
        },
    }

    if direction in ("citing", "both"):
        citing = _walk_citations(idx.biblio_db, work_id, "citing", depth)
        result["citing"] = citing
    if direction in ("cited_by", "both"):
        cited_by = _walk_citations(idx.biblio_db, work_id, "cited_by", depth)
        result["cited_by"] = cited_by

    return result


def _walk_citations(
    biblio: BiblioAuthority, work_id: str, direction: str, depth: int,
) -> List[Dict]:
    """BFS citation walk."""
    visited: set = set()
    frontier = [work_id]
    results: List[Dict] = []
    for d in range(depth):
        next_frontier: List[str] = []
        for wid in frontier:
            if wid in visited:
                continue
            visited.add(wid)
            rows = biblio.citing(wid) if direction == "citing" else biblio.cited_by(wid)
            for r in rows:
                r["depth"] = d + 1
                results.append(r)
                if r["work_id"] not in visited:
                    next_frontier.append(r["work_id"])
        frontier = next_frontier
    return results


@mcp.tool()
def resolve_reference(
    query: str,
    author: Optional[str] = None,
    year: Optional[int] = None,
) -> Dict:
    """Resolve a free-text bibliographic reference to a work in the
    authority database.

    Accepts forms like ``"Haeckel 1888"``, ``"Totton 1965 A Synopsis"``,
    or a raw citation string. Optionally pass ``author`` and ``year``
    separately for more precise matching. Returns the matched work with
    all known identifiers, in-corpus status, and citation counts.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return {"error": "bibliographic authority database not configured"}

    # Parse author and year from query if not provided separately
    if not author:
        import re
        # Try to extract "Author Year" or "Author, Year"
        m = re.match(r"^([A-Za-zÀ-ÿ\-\s]+?)[\s,]+(\d{4})\b", query.strip())
        if m:
            author = m.group(1).strip()
            if not year:
                year = int(m.group(2))

    if not author:
        return {"error": "could not parse author from query; pass author= explicitly"}

    # Title fragment is whatever remains after author+year
    title_frag = query
    if author:
        title_frag = title_frag.replace(author, "", 1).strip()
    if year:
        title_frag = title_frag.replace(str(year), "", 1).strip()
    title_frag = title_frag.strip(" ,;:")

    results = idx.biblio_db.search_works(
        author, year, title_frag if title_frag else None,
    )
    if not results:
        # Broaden: try without title
        results = idx.biblio_db.search_works(author, year)
    if not results:
        return {"not_found": True, "queried": query, "parsed_author": author, "parsed_year": year}
    if len(results) == 1:
        w = results[0]
        w["authors"] = idx.biblio_db.get_authors(w["work_id"])
        w["cited_by_count"] = idx.biblio_db.citation_count(w["work_id"])
        w["taxon_links"] = idx.biblio_db.taxon_links_for_work(w["work_id"])
        return w
    # Multiple matches — return all with citation counts
    for w in results:
        w["cited_by_count"] = idx.biblio_db.citation_count(w["work_id"])
    return {"matches": results, "count": len(results), "queried": query}


# #79 — warning footnote per provenance tier. Strings are public
# contract: the LLM is instructed to preserve them verbatim, and
# downstream document-quality tests grep for them. Don't change
# without coordinating both ends.
_PROVENANCE_WARNING = {
    "bib": "",
    "grobid_reconciled":
        "* generated via reconciliation in corpus, check if correct",
    "unresolved":
        "* reference not present in bibliography, check if correct",
}


@mcp.tool()
def format_citation(
    query: Optional[str] = None,
    work_id: Optional[str] = None,
    paper_hash: Optional[str] = None,
    style: str = "author-year",
) -> Dict[str, Any]:
    """Return a formatted citation string for a work in the corpus
    bibliographic authority DB (#79).

    **Use this for every citation you emit.** Never recombine
    author / year / journal / title in your own context — that's
    where amalgamated, hallucinated references come from. Paste
    ``formatted`` + ``inline`` verbatim; append non-empty ``warning``
    verbatim too.

    Resolve via exactly one of: ``query`` (free-text "Author Year
    [Title]"), ``work_id`` (DOI / corpus: / bhl:), or ``paper_hash``
    (12-hex SHA-256 prefix).

    Returns ``{work_id, formatted, inline, provenance, warning,
    fields}``. ``provenance`` ∈ {``bib`` (curated, no warning),
    ``grobid_reconciled`` (reconciliation footnote), ``unresolved``
    (bibliography-absence footnote)}.

    On ambiguous match returns ``{error: "ambiguous", matches: [...]}``;
    on no match returns ``{error: "not_found", ...}`` — say "not in
    the corpus" rather than fabricating one.
    """
    from bib.format import SUPPORTED_STYLES, format_citation as _format_str

    idx = _need_index()
    if idx.biblio_db is None:
        return {"error": "bibliographic authority database not configured"}
    if style not in SUPPORTED_STYLES:
        return {
            "error": f"unknown style {style!r}",
            "supported_styles": sorted(SUPPORTED_STYLES),
        }
    if sum(x is not None for x in (query, work_id, paper_hash)) != 1:
        return {
            "error": "provide exactly one of: query, work_id, paper_hash",
        }

    # Resolve to a single works row.
    work: Optional[Dict] = None
    if work_id is not None:
        work = idx.biblio_db.get_work(work_id)
        if work is None:
            return {"error": "not_found", "work_id": work_id}
    elif paper_hash is not None:
        work = idx.biblio_db.get_work_by_corpus_hash(paper_hash)
        if work is None:
            return {"error": "not_found", "paper_hash": paper_hash}
    else:  # query is not None
        import re
        m = re.match(r"^([A-Za-zÀ-ÿ\-\s]+?)[\s,]+(\d{4})\b", query.strip())
        if not m:
            return {
                "error": "could not parse author/year from query — "
                         "pass work_id explicitly or refine query",
                "queried": query,
            }
        author = m.group(1).strip()
        year = int(m.group(2))
        title_frag = query.replace(author, "", 1).strip()
        title_frag = title_frag.replace(str(year), "", 1).strip(" ,;:")
        results = idx.biblio_db.search_works(
            author, year, title_frag if title_frag else None,
        )
        if not results:
            # Broaden: drop the title constraint and retry.
            results = idx.biblio_db.search_works(author, year)
        if not results:
            return {
                "error": "not_found",
                "queried": query,
                "parsed_author": author,
                "parsed_year": year,
            }
        if len(results) > 1:
            return {
                "error": "ambiguous",
                "queried": query,
                "matches": [
                    {"work_id": r["work_id"], "title": r.get("title"),
                     "year": r.get("year"), "journal": r.get("journal")}
                    for r in results
                ],
            }
        work = idx.biblio_db.get_work(results[0]["work_id"])
        if work is None:  # search_works returned a row, get_work lost it
            return {"error": "not_found", "queried": query}

    # Assemble fields. volume + pages aren't on works.* — TODO when
    # the bibliography subsystem persists per-citation locator info.
    fields: Dict[str, Any] = {
        "authors": idx.biblio_db.get_authors(work["work_id"]),
        "year": work.get("year"),
        "title": work.get("title"),
        "journal": work.get("journal"),
        "volume": None,
        "pages": None,
        "doi": work.get("doi"),
    }
    rendered = _format_str(fields, style=style)
    provenance = idx.biblio_db.provenance(work["work_id"])
    return {
        "work_id": work["work_id"],
        "formatted": rendered["formatted"],
        "inline": rendered["inline"],
        "provenance": provenance,
        "warning": _PROVENANCE_WARNING[provenance],
        "fields": fields,
    }


@mcp.tool()
def get_missing_references(
    min_citations: int = 2,
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
    limit: int = 50,
) -> List[Dict]:
    """Works cited by corpus papers that are NOT in the corpus.

    Sorted by citation count (most-cited missing works first). Useful
    for identifying high-impact papers to add to the corpus. Filter by
    year range to focus on a particular era.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return [{"error": "bibliographic authority database not configured"}]

    query = """
        SELECT w.work_id, w.title, w.year, w.journal, w.doi,
               w.guid_type, COUNT(*) AS cited_by_count
        FROM citations c JOIN works w ON c.cited_work_id = w.work_id
        WHERE w.in_corpus = 0
    """
    params: list = []
    if year_from:
        query += " AND w.year >= ?"
        params.append(year_from)
    if year_to:
        query += " AND w.year <= ?"
        params.append(year_to)
    query += " GROUP BY c.cited_work_id HAVING COUNT(*) >= ?"
    params.append(int(min_citations))
    query += " ORDER BY cited_by_count DESC LIMIT ?"
    params.append(int(limit))

    cur = idx.biblio_db.conn.execute(query, params)
    results = []
    for row in cur:
        r = dict(row)
        r["authors"] = idx.biblio_db.get_authors(r["work_id"])
        results.append(r)
    return results


@mcp.tool()
def get_original_description(taxon_name: str) -> Dict:
    """Find the original description paper for a taxon.

    Resolves the taxon through DwC synonymy (``acceptedNameUsageID``),
    parses ``scientificNameAuthorship``, and looks up the corresponding
    work in the bibliographic authority database. Returns the work
    record, whether it's in the corpus, and if so the ``corpus_hash``
    for direct access via ``get_paper()``.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return {"error": "no taxonomy snapshot configured"}
    if idx.biblio_db is None:
        return {"error": "bibliographic authority database not configured"}

    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return {"not_found": True, "queried": taxon_name}

    aid = hit["accepted_taxon_id"]
    works = idx.biblio_db.work_for_taxon(aid)
    if not works:
        return {
            "taxon": hit,
            "original_description": None,
            "note": "no matching work found in the authority database",
        }

    # Enrich with authors and citation count
    for w in works:
        w["authors"] = idx.biblio_db.get_authors(w["work_id"])
        w["cited_by_count"] = idx.biblio_db.citation_count(w["work_id"])

    return {
        "taxon": hit,
        "original_description": works[0] if len(works) == 1 else None,
        "candidate_works": works if len(works) > 1 else None,
        "work": works[0] if len(works) == 1 else works,
    }


@mcp.tool()
def get_works_by_author(
    surname: str,
    in_corpus_only: bool = False,
    limit: int = 100,
) -> List[Dict]:
    """All works by an author across the full bibliographic authority
    database — not just corpus papers, but also cited references and
    stub works synthesized from taxonomic-authority strings.

    Each result includes ``cited_by_count`` (how many corpus papers
    cite it) and ``in_corpus`` flag. Use ``in_corpus_only=True`` to
    restrict to papers physically in the corpus.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return [{"error": "bibliographic authority database not configured"}]

    norm = surname.strip().lower()
    query = """
        SELECT DISTINCT w.work_id, w.title, w.year, w.journal, w.doi,
               w.in_corpus, w.corpus_hash, w.guid_type
        FROM works w JOIN work_authors wa ON w.work_id = wa.work_id
        WHERE wa.surname_normalized = ?
    """
    params: list = [norm]
    if in_corpus_only:
        query += " AND w.in_corpus = 1"
    query += " ORDER BY w.year LIMIT ?"
    params.append(int(limit))

    cur = idx.biblio_db.conn.execute(query, params)
    results = []
    for row in cur:
        r = dict(row)
        r["authors"] = idx.biblio_db.get_authors(r["work_id"])
        r["cited_by_count"] = idx.biblio_db.citation_count(r["work_id"])
        results.append(r)
    return results



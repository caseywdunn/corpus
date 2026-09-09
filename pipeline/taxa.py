"""Taxon + anatomy mention extraction over a Darwin Core taxonomy snapshot
and the anatomy lexicon.

See dev_docs/OVERVIEW.md "Taxonomic annotation" and "Lexicon-driven
annotation".

Two concerns:

1. :class:`TaxonomyDB` — thin read-only wrapper around the SQLite file
   produced by ``ingest_taxonomy.py``. Schema follows the Darwin Core
   Taxon class (``taxonID``, ``scientificName``, ``parentNameUsageID``,
   ``acceptedNameUsageID``, …) so any DwC source — a downloaded WoRMS,
   GBIF or iNaturalist export, or a user-supplied CSV — slots in behind
   the same interface.
2. :func:`extract_taxon_mentions` and :func:`extract_lexicon_mentions`
   scan chunks for candidate name spans and lexicon terms respectively,
   returning flat lists of mentions with chunk_id + char offsets and
   per-taxon / per-term rollups.

Why not gnfinder: a simple regex-based candidate extractor + name lookup
covers the 80% case (explicit binomials, trinomials, and uninomials)
without adding a Go binary or a network-dependent name-verification API.
When we find historical-spelling misses, adding a gnfinder stage is a
clean extension — the interface here is name → taxonomy resolution, not
gnfinder-specific.
"""

from __future__ import annotations

import logging
import re
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Taxonomy SQLite access (Darwin Core Taxon schema)
# ---------------------------------------------------------------------------


class TaxonomyDB:
    """Read-only lookup into a Darwin Core taxonomy SQLite snapshot.

    Schema is produced by ``ingest_taxonomy.py`` (any source: a DwC file,
    a Darwin Core Archive, or the WoRMS REST API). Field names match DwC
    Taxon-class terms verbatim:

    * ``taxon_id``                       — DwC ``taxonID`` (TEXT)
    * ``scientific_name``                — DwC ``scientificName``
    * ``scientific_name_authorship``     — DwC ``scientificNameAuthorship``
    * ``taxon_rank``                     — DwC ``taxonRank``
    * ``taxonomic_status``               — DwC ``taxonomicStatus``
    * ``parent_name_usage_id``           — DwC ``parentNameUsageID``
    * ``accepted_name_usage_id``         — DwC ``acceptedNameUsageID``

    A single DB connection is held; the class is safe to use from one
    process. For parallel workers, each worker should open its own
    TaxonomyDB — sqlite3 connections aren't thread-safe by default.

    Methods are designed around the two queries the pipeline needs:

    * ``lookup(name)`` — given a candidate text span, resolve to the
      accepted taxon (following ``acceptedNameUsageID`` synonymy).
      Returns None if the name isn't in the snapshot.
    * ``name_set()`` — the set of lowercased names in the snapshot. The
      caller uses this to pre-filter candidate spans before calling
      ``lookup()``, so we don't hit SQLite for every capitalized word in
      the corpus.
    """

    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(
                f"Taxonomy SQLite not found at {self.db_path}. "
                f"Build it with: python -m pipeline.taxonomy_ingest --source <dwc|dwca|worms> ..."
            )
        # uri=... read-only for safety; no writer process should touch
        # the snapshot during pipeline runs.
        self.conn = sqlite3.connect(
            f"file:{self.db_path}?mode=ro", uri=True, check_same_thread=False
        )
        self.conn.row_factory = sqlite3.Row
        self._name_set_cache: Optional[Set[str]] = None

    def close(self) -> None:
        self.conn.close()

    def __enter__(self) -> "TaxonomyDB":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # ---- lookup ----

    def name_set(self) -> Set[str]:
        """Return the full lowercased name set. Cached on first call."""
        if self._name_set_cache is None:
            cur = self.conn.execute("SELECT DISTINCT name_lowercase FROM names")
            self._name_set_cache = {row[0] for row in cur}
        return self._name_set_cache

    def lookup(self, name: str) -> Optional[Dict]:
        """Resolve a text span to an accepted DwC taxon.

        Returns a dict with the matched-name's taxon_id and the accepted
        taxon's ``taxon_id`` / ``scientific_name`` (following
        ``accepted_name_usage_id``; may be the same record), plus
        authorship / rank / status. Returns None if not found.

        Lookup is case-insensitive. When multiple taxa share a lowercased
        name, prefer an accepted primary name, then an unaccepted primary
        name, then a synonym alias; ties are ordered by taxon ID.  This keeps
        the result independent of SQLite insertion order while preserving a
        directly named unresolved taxon rather than silently forcing it onto
        a homonymous synonym target.
        """
        key = (name or "").strip().lower()
        if not key:
            return None
        cur = self.conn.execute(
            """
            SELECT n.name, n.name_type, n.taxon_id,
                   t.scientific_name AS matched_name,
                   t.scientific_name_authorship AS authorship,
                   t.taxon_rank AS rank,
                   t.taxonomic_status AS status,
                   t.accepted_name_usage_id AS accepted_id,
                   t.accepted_name
            FROM names n
            JOIN taxa t ON t.taxon_id = n.taxon_id
            WHERE n.name_lowercase = ?
            ORDER BY CASE n.name_type
                       WHEN 'accepted' THEN 0
                       WHEN 'unaccepted' THEN 1
                       ELSE 2
                     END,
                     n.taxon_id,
                     n.name
            """,
            (key,),
        )
        rows = list(cur)
        if not rows:
            return None

        row = rows[0]
        accepted_id = row["accepted_id"] or row["taxon_id"]
        accepted_name = row["accepted_name"] or row["matched_name"]
        # If the matched taxon points at an accepted_id we don't have
        # locally (rare — pruning should keep both ends of a synonym
        # link), fall back to the matched record's own scientific name.
        if accepted_id != row["taxon_id"]:
            cur2 = self.conn.execute(
                "SELECT scientific_name, scientific_name_authorship, taxon_rank "
                "FROM taxa WHERE taxon_id = ?",
                (accepted_id,),
            )
            acc = cur2.fetchone()
            if acc:
                accepted_name = acc["scientific_name"] or accepted_name
        return {
            "matched_taxon_id": row["taxon_id"],
            "matched_name": row["matched_name"],
            "name_type": row["name_type"],
            "accepted_taxon_id": accepted_id,
            "accepted_name": accepted_name,
            "authorship": row["authorship"],
            "rank": row["rank"],
            "status": row["status"],
        }


# ---------------------------------------------------------------------------
# Taxon-name candidate extraction from text
# ---------------------------------------------------------------------------

# Match sequences of Latin-letter tokens starting with a capital, with up to
# two lowercase-initial tokens after (to cover trinomials like
# "Bargmannia elongata maculata"). The opening `\b` ensures we don't match
# inside a word; `(?:-|\s)` allows hyphenated species names. Excludes
# leading non-letter by using \A-less \b anchor.
_NAME_CANDIDATE_RE = re.compile(
    r"""
    \b
    (
      [A-Z][a-z]{2,}                    # Genus — at least 3 chars to skip initials
      (?:
        (?:\s|-)
        [a-z][a-z\-]{2,}                # species epithet
        (?:
          (?:\s|-)
          [a-z][a-z\-]{2,}              # subspecies
        )?
      )?
    )
    \b
    """,
    re.VERBOSE,
)


# Common words we should never treat as uninomial taxon candidates even if
# they happen to match a taxonomy name (usually at a higher rank). These
# are English words that collide with marine taxon names and cause obvious
# false positives in biological prose.
_COMMON_WORD_STOPLIST = frozenset(
    s.lower() for s in [
        "Figure", "Figures", "Fig", "Table", "Tables", "Plate", "Plates",
        "Section", "Sections", "Chapter", "Chapters",
        "Abstract", "Introduction", "Methods", "Results", "Discussion",
        "Acknowledgements", "References", "Appendix",
        "January", "February", "March", "April", "May", "June",
        "July", "August", "September", "October", "November", "December",
        "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday",
        "University", "Museum", "Institute", "Department", "Laboratory",
        "Pacific", "Atlantic", "Indian", "Arctic", "Antarctic",  # used generically
        "North", "South", "East", "West",
        "Marine", "Ocean", "Sea", "Bay", "Gulf",
        "Specimen", "Specimens", "Sample", "Samples",
        "Species", "Genus", "Family", "Order", "Class",
    ]
)


def _extract_name_candidates(text: str) -> Iterable[Tuple[str, int, int]]:
    """Yield (name, start, end) tuples for capitalized candidate spans."""
    for m in _NAME_CANDIDATE_RE.finditer(text):
        name = m.group(1)
        if name.lower() in _COMMON_WORD_STOPLIST:
            continue
        yield name, m.start(1), m.end(1)


# An abbreviated genus followed by an epithet: `Ph. pelagica`, `P. physalis`,
# and — because this corpus is OCR'd — `Ph, pelagica` with a comma for the
# period. The prefix is 1-4 letters; the epithet is spelled out, which is
# what makes the expansion checkable at all.
#
# This regex is deliberately loose, because it is not what decides
# anything: `Taf. iii` and `No. species` match it too. Two gates behind it
# do the deciding — the prefix must extend a genus written out in full
# elsewhere in the *same document*, and the expansion must resolve in the
# taxonomy snapshot. A bibliographic abbreviation clears neither.
_ABBREV_BINOMIAL_RE = re.compile(
    r"""
    \b
    ([A-Z][a-z]{0,3})           # abbreviated genus: Ph, P, Phys
    \s* [.,] \s*                # period, or an OCR'd comma
    ([a-z][a-z\-]{2,})          # epithet, spelled out
    \b
    """,
    re.VERBOSE,
)


def _genus_of(name: str) -> str:
    """The genus token of a resolved name span."""
    return re.split(r"[\s-]", name.strip(), 1)[0]


def _expand_abbreviation(
    prefix: str, epithet: str, genera_in_full: "Counter[str]",
    taxonomy: "TaxonomyDB", name_set: Set[str],
) -> Tuple[Optional[str], Optional[Dict], List[str]]:
    """Resolve ``prefix. epithet`` against one document's own genera.

    Returns ``(expanded_name, resolved, candidates)``. The first two are
    set only when the candidates agree on a single taxon; ``candidates``
    is every binomial that cleared both gates, so a genuinely ambiguous
    case can be reported rather than guessed at.

    Three gates. The genus must be one this document writes out in full,
    the expansion must be a name in the taxonomy snapshot, and the
    survivors must resolve to one accepted taxon.

    The second gate is what makes this safe rather than a guess. `Ph.` is
    genuinely ambiguous in a siphonophore corpus — *Physalia* and
    *Physophora* are both in it — so document context alone would be
    guessing. But the taxonomy knows `Physalia pelagica` is a name and
    `Physophora pelagica` is not, and the epithet is printed right there.

    The third gate exists because the second over-reports. Historical
    spellings live in the snapshot as names of their own: Olfers 1824
    prints both *Physalia* and *Physalis*, and `Physalia pelagica` and
    `Physalis pelagica` are two names for accepted taxon 135479. Counting
    name strings called that ambiguous and dropped a mention that was
    never in doubt. Ambiguity is disagreement about the *taxon*.

    Where several spellings do agree, the document's own preference picks
    the representative — Olfers writes *Physalia* twice and *Physalis*
    once — so the recorded name is the one the author mostly used.
    """
    lowered = prefix.lower()
    candidates = [
        f"{genus} {epithet}"
        for genus, _ in sorted(genera_in_full.most_common(),
                               key=lambda kv: (-kv[1], kv[0]))
        if genus.lower().startswith(lowered)
        and f"{genus} {epithet}".lower() in name_set
    ]
    if not candidates:
        return None, None, []
    resolved_by_taxon: Dict[str, Tuple[str, Dict]] = {}
    for name in candidates:
        r = taxonomy.lookup(name)
        if r is None:
            continue
        resolved_by_taxon.setdefault(r["accepted_taxon_id"], (name, r))
    if len(resolved_by_taxon) != 1:
        return None, None, candidates
    name, resolved = next(iter(resolved_by_taxon.values()))
    return name, resolved, candidates


# taxa.json ships in the served bundle, so the ambiguity report is a
# bounded diagnostic rather than a parallel mention list.
_MAX_UNRESOLVED_RECORDED = 50


def extract_taxon_mentions(
    chunks: List[Dict],
    taxonomy: TaxonomyDB,
    *,
    min_chars: int = 4,
) -> Dict:
    """Scan each chunk for taxon names present in the taxonomy snapshot.

    Returns a dict with ``mentions`` (flat, ordered by chunk then
    position) and ``taxa`` (rolled up one row per accepted ``taxon_id``
    with its mention count). Both forms are useful: the flat list for
    per-chunk or locality-resolving queries; the rollup for corpus-wide
    "which species appear in which papers".

    Only mentions whose ``accepted_taxon_id`` resolves to a taxon in the
    snapshot are recorded — this filters out generic words that happen
    to collide with taxonomic names outside the configured subtree.

    Runs in two passes over the document (#164). The first resolves names
    written out in full and, as a side effect, learns which genera this
    document spells out. The second expands abbreviated binomials —
    `Ph. pelagica` — against that set. Taxonomic literature abbreviates
    the genus after first mention, so for a corpus of original
    descriptions this is the central case rather than an edge one: the
    paper that *erects* a species is the one least likely to spell the
    genus out on every line. Olfers 1824 is a five-species key for
    *Physalia* that yielded one genus-level taxon and no species.

    Expansions are recorded with ``method="abbreviated_genus"`` and keep
    the printed form in ``mention_text``, so nothing downstream has to
    take an inferred name for an observed one. Abbreviations the pass
    declined to resolve land in ``abbreviations_unresolved`` rather than
    disappearing.
    """
    name_set = taxonomy.name_set()
    if not name_set:
        logger.warning("Taxonomy name set is empty; no taxon mentions will be recorded")

    # Per-chunk so pass 2 can interleave its mentions in text order
    # rather than appending a block at the end.
    by_chunk: List[List[Dict]] = []
    genera_in_full: "Counter[str]" = Counter()

    for ch in chunks:
        chunk_mentions: List[Dict] = []
        by_chunk.append(chunk_mentions)
        text = ch.get("text", "") or ""
        if not text:
            continue
        for name, start, end in _extract_name_candidates(text):
            if len(name) < min_chars:
                continue
            # Try progressively shorter prefixes: full trinomial, then
            # binomial, then genus. The regex can greedily swallow
            # trailing lowercase words ("Agalma elegans and"), so the
            # full match often doesn't hit name_set — but a prefix will.
            # Longest match wins so that "Bargmannia elongata" is
            # recorded as the species, not just the genus.
            tokens = re.split(r"(\s|-)", name)  # preserve separators
            words = [t for i, t in enumerate(tokens) if i % 2 == 0]
            word_count = len(words)

            resolved = None
            matched_text = None
            for n_words in range(word_count, 0, -1):
                take = n_words * 2 - 1  # words + separators between them
                candidate_text = "".join(tokens[:take])
                if candidate_text.lower() in name_set:
                    resolved = taxonomy.lookup(candidate_text)
                    matched_text = candidate_text
                    end = start + len(candidate_text)
                    break

            if resolved is None:
                continue

            genera_in_full[_genus_of(matched_text)] += 1
            chunk_mentions.append(
                {
                    "chunk_id": ch.get("chunk_id"),
                    "text_span": [start, end],
                    "matched_text": matched_text,
                    "matched_taxon_id": resolved["matched_taxon_id"],
                    "name_type": resolved["name_type"],
                    "accepted_taxon_id": resolved["accepted_taxon_id"],
                    "accepted_name": resolved["accepted_name"],
                    "authorship": resolved["authorship"],
                    "rank": resolved["rank"],
                }
            )

    # --- Pass 2: abbreviated binomials, against this document's genera ---
    unresolved: List[Dict] = []
    n_expanded = 0
    for ch, chunk_mentions in zip(chunks, by_chunk):
        text = ch.get("text", "") or ""
        if not text or not genera_in_full:
            continue
        # Spans pass 1 already claimed. An abbreviation cannot overlap a
        # name written out in full, and `P. Sars` style author initials
        # next to a resolved name must not be re-read as an epithet.
        claimed = [tuple(m["text_span"]) for m in chunk_mentions]
        for m in _ABBREV_BINOMIAL_RE.finditer(text):
            start, end = m.start(0), m.end(0)
            if any(s < end and start < e for s, e in claimed):
                continue
            prefix, epithet = m.group(1), m.group(2)
            expanded, resolved, candidates = _expand_abbreviation(
                prefix, epithet, genera_in_full, taxonomy, name_set,
            )
            if expanded is None or resolved is None:
                if candidates:
                    unresolved.append({
                        "chunk_id": ch.get("chunk_id"),
                        "text_span": [start, end],
                        "mention_text": m.group(0),
                        "candidates": candidates,
                        "reason": "ambiguous_abbreviation",
                    })
                continue
            n_expanded += 1
            chunk_mentions.append({
                "chunk_id": ch.get("chunk_id"),
                "text_span": [start, end],
                "matched_text": expanded,
                # What is actually printed on the page. The schema has
                # kept these apart all along; only the writer collapsed
                # them.
                "mention_text": m.group(0),
                "matched_taxon_id": resolved["matched_taxon_id"],
                "name_type": resolved["name_type"],
                "accepted_taxon_id": resolved["accepted_taxon_id"],
                "accepted_name": resolved["accepted_name"],
                "authorship": resolved["authorship"],
                "rank": resolved["rank"],
                "method": "abbreviated_genus",
                "expanded_from": prefix,
            })

    mentions: List[Dict] = []
    for chunk_mentions in by_chunk:
        mentions.extend(sorted(chunk_mentions, key=lambda m: m["text_span"][0]))

    taxa_rollup: Dict[str, Dict] = {}
    for m in mentions:
        accepted_id = m["accepted_taxon_id"]
        bucket = taxa_rollup.setdefault(
            accepted_id,
            {
                "accepted_taxon_id": accepted_id,
                "accepted_name": m["accepted_name"],
                "authorship": m["authorship"],
                "rank": m["rank"],
                "mention_count": 0,
                "first_chunk": m["chunk_id"],
            },
        )
        bucket["mention_count"] += 1

    taxa_list = sorted(
        taxa_rollup.values(),
        key=lambda r: (-r["mention_count"], r["accepted_name"] or ""),
    )
    out = {
        "total_mentions": len(mentions),
        "unique_taxa": len(taxa_list),
        "mentions": mentions,
        "taxa": taxa_list,
        "abbreviations_expanded": n_expanded,
        # Bounded: this is a diagnostic, not a second mention list, and
        # taxa.json ships in the served bundle.
        "abbreviations_unresolved_count": len(unresolved),
        "abbreviations_unresolved": unresolved[:_MAX_UNRESOLVED_RECORDED],
    }
    if len(unresolved) > _MAX_UNRESOLVED_RECORDED:
        logger.info(
            "%d ambiguous genus abbreviations, recording the first %d",
            len(unresolved), _MAX_UNRESOLVED_RECORDED,
        )
    return out


# ---------------------------------------------------------------------------
# Anatomy lexicon and term extraction
# ---------------------------------------------------------------------------


RESERVED_ARTIFACT_STEMS = frozenset({
    "taxa", "summary", "metadata", "references", "text", "chunks", "figures",
    "intext_citations", "pipeline_state", "scan_detection", "docling_doc",
    "annotation_outputs", "figure_materialization_base", "grobid.tei.provenance",
})


def validate_category(category):
    """A category may name its own file, never a core artifact or a path."""
    if not category or Path(category).name != category or "\\" in category or category in {".", ".."} or category in RESERVED_ARTIFACT_STEMS:
        raise ValueError(f"Invalid or reserved lexicon category: {category!r}")


def load_lexicon(path: Path) -> Dict[str, Dict[str, Dict]]:
    """Load a multi-category lexicon YAML.

    The file is two-level: each top-level key is a category name
    (``anatomy``, ``biogeography``, …) whose value is the per-term
    mapping that drives extraction:

        anatomy:
          pneumatophore:
            synonyms: [pneumatophores, float]
            translations: {de: [Luftblase]}
            description: Apical gas-filled float.
          nectosome:
            ...
        biogeography:
          pelagic:
            synonyms: [open water]
            description: ...

    Returns ``{category: {canonical_term: {synonyms, translations,
    description}}}``. Category names are normalized to lowercase. Term
    keys are taken as-written (lowercase snake_case recommended).
    """
    import yaml
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(
            f"{path}: lexicon root must be a mapping of category → terms"
        )
    out: Dict[str, Dict[str, Dict]] = {}
    for category, terms in data.items():
        category = str(category).strip().lower()
        if not category:
            continue
        validate_category(category)
        if not isinstance(terms, dict):
            raise ValueError(
                f"{path}: category '{category}' must be a mapping of "
                f"term → {{synonyms, translations, description}}"
            )
        section: Dict[str, Dict] = {}
        for canonical, entry in terms.items():
            entry = entry or {}
            if not isinstance(entry, dict):
                raise ValueError(
                    f"{path}: category '{category}', term '{canonical}': "
                    f"value must be a mapping of {{synonyms, translations, "
                    f"description}} (got {type(entry).__name__}). If you "
                    f"meant a list of synonyms, write "
                    f"`{canonical}: {{synonyms: [...]}}`."
                )
            section[canonical] = {
                "synonyms": list(entry.get("synonyms") or []),
                "translations": dict(entry.get("translations") or {}),
                "description": entry.get("description") or "",
            }
        out[category] = section
    return out


def lexicon_fingerprints(path: Path) -> Dict[str, Dict[str, object]]:
    """Per-category content hashes for the lexicon at ``path``.

    Returns ``{category: {"path": ..., "sha256": ..., "size": ...}}``.
    The SHA is taken over the *canonical* JSON of just that category's
    section, so a change to one category does not perturb another's
    fingerprint. ``path`` is recorded so corpus_status can attribute
    a stale fingerprint back to its source file.
    """
    import hashlib
    import json as _json
    import yaml
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    out: Dict[str, Dict[str, object]] = {}
    for category, section in (data or {}).items():
        category = str(category).strip().lower()
        if not category:
            continue
        canonical = _json.dumps(
            section or {}, sort_keys=True, ensure_ascii=False,
        ).encode("utf-8")
        out[category] = {
            "path": str(path),
            "sha256": hashlib.sha256(canonical).hexdigest(),
            "size": len(canonical),
        }
    return out


# Noun-inflection endings applied to `translations` forms, per language
# (#165). A curated list rather than a general stemmer: every generated
# form is added to the variant map explicitly, so matching stays
# whole-word exact and the map stays greppable.
#
# `de`, `fr` and `ru` are measured — the endings below are the ones that
# actually follow a lexicon stem in the reference library, with the
# per-ending hit counts that justified each. The others carry that
# language's ordinary noun plurals and are there so a lexicon extended
# to them is not silently worse off; they have no corpus evidence yet.
#
# There is deliberately no `en` entry, and it is not an omission. English
# variants are hand-listed in `synonyms`, which is what the lexicon's own
# documentation tells curators to do, and the survey shows why: suffixing
# English stems matches `Cnidaria` 4,444 times and `cnidarian(s)` 3,740
# more from `cnida`, which is the phylum rather than the nematocyst, plus
# `stemmed` from `stem`, `floating` from `float` and `siphoning` from
# `siphon`. Those are wrong, not merely loose.
_INFLECTION_SUFFIXES: Dict[str, Tuple[str, ...]] = {
    # Schwimmglocken 4,623 · Deckstücke 663 · Deckstücken 242 ·
    # Deckstückes 225 · Deckstücks/Tentakels 988
    "de": ("n", "en", "e", "es", "s", "ns"),
    # tentacules/cnidocytes/nématocystes/bractées 3,168 · palpones 6
    "fr": ("s", "es", "x"),
    # нектофора 831 · нектофоры 603 · нектофоров 406 · нектофорами 77 ·
    # пневматофором 70 · нектофорам 18 · нектофорах 16 · нектофору 14 ·
    # нектофоре 13
    "ru": ("а", "я", "ы", "и", "ов", "ев", "ей", "ам", "ям", "ами", "ями",
           "ах", "ях", "ом", "ем", "у", "ю", "е", "ой"),
    "es": ("s", "es"),
    "pt": ("s", "es"),
    "it": ("i", "e"),
    "nl": ("n", "en", "s"),
    "la": ("e", "ae", "a", "um", "i", "is", "es", "orum", "arum"),
}

# Languages whose case endings replace a stem-final vowel rather than
# stacking on it. Russian `личинка` declines to `личинки`, not
# `личинкаи`, so the ending has to be appended to `личинк` as well.
# Appending to the bare stem still covers the consonant-final majority
# (`нектофор` → `нектофора`), so both forms are generated.
#
# What suffixing cannot reach, at any ending-list length, is an
# inflection that changes the stem: Russian genitive plurals insert a
# fill vowel (`личинка` → `личин|о|к`) and German umlaut plurals change
# one (`Saugmagen` → `Saugmägen`, `Fangfaden` → `Fangfäden`, both printed
# by Eschscholtz). Those want a per-language morphology table, or the
# curator listing the form under `synonyms`, which works today.
# `tests/test_lexicon_inflection.py` pins the gap so it stays recorded
# rather than merely absent.
_VOWEL_DROPPING = {"ru"}
_RU_VOWELS = "аеёиоуыэюя"

# Below this many characters a stem plus an ending is as likely to be an
# unrelated short word as an inflection. The shortest real translation
# form in the reference lexicon is `Larve` at 5.
_MIN_INFLECTABLE_STEM = 5


def _inflected_forms(form: str, lang: str) -> List[str]:
    """Surface forms of ``form`` a reader of ``lang`` might have printed.

    Returns the empty list for languages with no ending table and for
    stems too short to suffix safely.
    """
    suffixes = _INFLECTION_SUFFIXES.get((lang or "").lower())
    if not suffixes or len(form) < _MIN_INFLECTABLE_STEM:
        return []
    # A multi-word translation inflects on its head, not by having an
    # ending stuck on the end of the phrase, so leave those alone.
    if " " in form.strip() or "-" in form:
        return []
    stems = [form]
    if (lang or "").lower() in _VOWEL_DROPPING and form[-1].lower() in _RU_VOWELS:
        stems.append(form[:-1])
    out: List[str] = []
    for stem in stems:
        for suffix in suffixes:
            out.append(stem + suffix)
    return out


def _build_lexicon_matcher(lexicon: Dict[str, Dict]) -> Tuple[re.Pattern, Dict[str, str]]:
    """Compile one big alternation regex and a variant→canonical map.

    Variants are matched as whole words, case-insensitive. Canonical
    names (e.g., ``nectophore``) ARE matched — their own keys count as
    variants of themselves so you don't have to repeat them in the
    ``synonyms`` list.

    Non-English translations are additionally expanded through
    :func:`_inflected_forms`, because an enumerated surface-form set is
    the wrong shape for an inflecting language: Eschscholtz prints
    `Luftblasen` where the lexicon lists `Luftblase`, and Vanhöffen 1906
    prints `Schwimmglocken` 41 times against 7 of `Schwimmglocke` — so
    the base-form-only match found a minority of its own mentions and a
    German paper could report anatomy coverage of exactly zero, which
    reads as "nothing here" rather than "not indexed" (#165).

    An explicitly curated form always wins over a generated one, so a
    lexicon can correct a bad generation by listing the right form.
    """
    variant_to_canonical: Dict[str, str] = {}
    generated: Dict[str, str] = {}
    for canonical, entry in lexicon.items():
        canonical_display = canonical.replace("_", " ")
        variants = {canonical, canonical_display, *entry["synonyms"]}
        for lang, lang_variants in entry.get("translations", {}).items():
            variants.update(lang_variants)
            for form in lang_variants or []:
                for inflected in _inflected_forms(form, lang):
                    generated.setdefault(inflected.lower(), canonical)
        for v in variants:
            if v:
                variant_to_canonical[v.lower()] = canonical
    # Curated forms are already in the map; only fill the gaps, so a
    # generated form can never displace a term someone typed on purpose.
    for variant, canonical in generated.items():
        variant_to_canonical.setdefault(variant, canonical)

    # Sort longest-first so multi-word variants match before shorter ones
    # (e.g., "swimming bell" before "bell").
    escaped = sorted(
        (re.escape(v) for v in variant_to_canonical),
        key=len,
        reverse=True,
    )
    pattern = re.compile(
        r"(?<![\w])(" + "|".join(escaped) + r")(?![\w])",
        re.IGNORECASE,
    )
    return pattern, variant_to_canonical


def extract_lexicon_mentions(
    chunks: List[Dict],
    lexicon: Dict[str, Dict],
    *,
    category: Optional[str] = None,
) -> Dict:
    """Scan chunks for a domain lexicon's terms.

    Generalization of the original anatomy-only extractor (#24): the
    same matcher produces mentions for any user-curated category
    (anatomy, biogeography, life history, …). When ``category`` is
    supplied it lands on the output dict so consumers can tell
    multiple categories apart.

    Same output shape as :func:`extract_taxon_mentions` — flat mentions
    list plus a per-canonical-term rollup. Matches are case-insensitive
    and whole-word; see :func:`_build_lexicon_matcher`.
    """
    if not lexicon:
        out_empty = {"total_mentions": 0, "unique_terms": 0, "mentions": [], "terms": []}
        if category is not None:
            out_empty["category"] = category
        return out_empty

    pattern, variant_to_canonical = _build_lexicon_matcher(lexicon)

    mentions: List[Dict] = []
    term_counts: Counter = Counter()
    term_first_chunk: Dict[str, str] = {}

    for ch in chunks:
        text = ch.get("text", "") or ""
        if not text:
            continue
        for m in pattern.finditer(text):
            variant = m.group(1)
            canonical = variant_to_canonical.get(variant.lower())
            if canonical is None:
                continue
            mentions.append(
                {
                    "chunk_id": ch.get("chunk_id"),
                    "text_span": [m.start(1), m.end(1)],
                    "matched_text": variant,
                    "canonical": canonical,
                }
            )
            term_counts[canonical] += 1
            term_first_chunk.setdefault(canonical, ch.get("chunk_id"))

    terms_list = [
        {
            "canonical": canonical,
            "mention_count": count,
            "first_chunk": term_first_chunk.get(canonical),
            "description": lexicon.get(canonical, {}).get("description", ""),
        }
        for canonical, count in term_counts.most_common()
    ]
    out = {
        "total_mentions": len(mentions),
        "unique_terms": len(terms_list),
        "mentions": mentions,
        "terms": terms_list,
    }
    if category is not None:
        out["category"] = category
    return out

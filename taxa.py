"""Taxon + anatomy mention extraction over the WoRMS snapshot and the
anatomy lexicon.

Phase F. See PLAN.md §3 ("Corpus-level inputs") and §4 ("Taxonomic
mention extraction + resolution").

Two concerns:

1. :class:`WormsDB` — thin read-only wrapper around the SQLite file
   produced by ``scripts/ingest_worms.py``. Case-insensitive name lookup
   via the indexed ``names.name_lowercase`` column resolves to an
   AphiaID + the accepted-name AphiaID (following synonymy).
2. :func:`extract_taxon_mentions` and :func:`extract_anatomy_mentions`
   scan chunks for candidate name spans and anatomy terms respectively,
   returning flat lists of mentions with chunk_id + char offsets and
   per-taxon / per-term rollups.

Why not gnfinder: for the initial Phase F we use a simple regex-based
candidate extractor + WoRMS name lookup. This covers the 80% case
(explicit binomials, trinomials, and uninomials) without adding a Go
binary or a network-dependent name-verification API. When we find
historical-spelling misses on the Pugh-Phil library, adding a gnfinder
stage is a clean extension — the interface here is name→WoRMS resolution,
not gnfinder-specific.
"""

from __future__ import annotations

import logging
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# WoRMS SQLite access
# ---------------------------------------------------------------------------


class WormsDB:
    """Read-only lookup into the WoRMS SQLite snapshot.

    A single DB connection is held; the class is safe to use from one
    process (the per-PDF loop is serial). For parallel workers, each
    worker should open its own WormsDB — sqlite3 connections aren't
    thread-safe by default.

    Methods are designed around the two queries the pipeline needs:

    * ``lookup(name)`` — given a candidate text span, resolve to the
      accepted WoRMS taxon (following synonymy). Returns None if the
      name isn't in the snapshot.
    * ``name_set()`` — the set of lowercased names in the snapshot. The
      caller uses this to pre-filter candidate spans before calling
      lookup(), so we don't hit SQLite for every capitalized word in the
      corpus.
    """

    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(
                f"WoRMS SQLite not found at {self.db_path}. "
                f"Run: python scripts/ingest_worms.py"
            )
        # uri=... read-only for safety; no writer process should ever touch
        # the snapshot during pipeline runs.
        self.conn = sqlite3.connect(
            f"file:{self.db_path}?mode=ro", uri=True, check_same_thread=False
        )
        self.conn.row_factory = sqlite3.Row
        self._name_set_cache: Optional[Set[str]] = None

    def close(self) -> None:
        self.conn.close()

    def __enter__(self) -> "WormsDB":
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
        """Resolve a text span to an accepted WoRMS taxon.

        Returns a dict with AphiaIDs and canonical names for both the
        matched name and the accepted name it points to (may be the
        same), plus rank / authority / status. Returns None if not found.

        The lookup is case-insensitive. When multiple taxa share a
        lowercased name (rare but possible — homonyms across kingdoms),
        the first accepted match is preferred, then the first match
        outright.
        """
        key = (name or "").strip().lower()
        if not key:
            return None
        cur = self.conn.execute(
            """
            SELECT n.name, n.name_type, n.aphia_id,
                   t.scientific_name AS matched_name,
                   t.authority, t.rank, t.status,
                   t.valid_aphia_id, t.valid_name
            FROM names n
            JOIN taxa t ON t.aphia_id = n.aphia_id
            WHERE n.name_lowercase = ?
            """,
            (key,),
        )
        rows = list(cur)
        if not rows:
            return None

        # Prefer an 'accepted' match when there are multiple.
        row = next(
            (r for r in rows if r["name_type"] == "accepted"),
            rows[0],
        )
        accepted_id = row["valid_aphia_id"] or row["aphia_id"]
        accepted_name = row["valid_name"] or row["matched_name"]
        # If the matched taxon points at a valid_aphia_id we don't have
        # locally (should be rare — the walk covers all descendants), fall
        # back to the matched record's own name.
        if accepted_id != row["aphia_id"]:
            cur2 = self.conn.execute(
                "SELECT scientific_name, authority, rank FROM taxa WHERE aphia_id = ?",
                (accepted_id,),
            )
            acc = cur2.fetchone()
            if acc:
                accepted_name = acc["scientific_name"] or accepted_name
        return {
            "matched_aphia_id": row["aphia_id"],
            "matched_name": row["matched_name"],
            "name_type": row["name_type"],
            "accepted_aphia_id": accepted_id,
            "accepted_name": accepted_name,
            "authority": row["authority"],
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
# they happen to match a WoRMS name (usually at a higher rank). These are
# English words that collide with marine taxon names and cause obvious
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
        "Figure",  # repeat intentional so extending is cheap
    ]
)


def _extract_name_candidates(text: str) -> Iterable[Tuple[str, int, int]]:
    """Yield (name, start, end) tuples for capitalized candidate spans."""
    for m in _NAME_CANDIDATE_RE.finditer(text):
        name = m.group(1)
        if name.lower() in _COMMON_WORD_STOPLIST:
            continue
        # Strip trailing punctuation that wasn't consumed by `\b` — e.g.,
        # "Marrus claudanielis." shouldn't include the period.
        yield name, m.start(1), m.end(1)


def extract_taxon_mentions(
    chunks: List[Dict],
    worms: WormsDB,
    *,
    min_chars: int = 4,
) -> Dict:
    """Scan each chunk for taxon names present in the WoRMS snapshot.

    Returns a dict with ``mentions`` (flat, ordered by chunk then
    position) and ``taxa`` (rolled up one row per accepted AphiaID with
    its mention count). Both forms are useful: the flat list for
    per-chunk or locality-resolving queries; the rollup for corpus-wide
    "which species appear in which papers".

    Only mentions whose ``accepted_aphia_id`` resolves to a taxon in our
    snapshot are recorded — this filters out generic words that happen
    to collide with WoRMS names outside the Siphonophorae subtree.
    """
    name_set = worms.name_set()
    if not name_set:
        logger.warning("WoRMS name set is empty; no taxon mentions will be recorded")

    mentions: List[Dict] = []
    taxa_rollup: Dict[int, Dict] = {}

    for ch in chunks:
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
            # Rebuild candidate names from the longest to shortest word
            # prefix. With preserved separators, tokens look like
            # ["Agalma", " ", "elegans", " ", "and"]; we slice by word
            # index. Words are at even positions 0, 2, 4.
            words = [t for i, t in enumerate(tokens) if i % 2 == 0]
            word_count = len(words)

            resolved = None
            matched_text = None
            for n_words in range(word_count, 0, -1):
                # Reassemble the first n_words words with their separators.
                take = n_words * 2 - 1  # words + separators between them
                candidate_text = "".join(tokens[:take])
                if candidate_text.lower() in name_set:
                    resolved = worms.lookup(candidate_text)
                    matched_text = candidate_text
                    end = start + len(candidate_text)
                    break

            if resolved is None:
                continue

            mentions.append(
                {
                    "chunk_id": ch.get("chunk_id"),
                    "text_span": [start, end],
                    "matched_text": matched_text,
                    "matched_aphia_id": resolved["matched_aphia_id"],
                    "name_type": resolved["name_type"],
                    "accepted_aphia_id": resolved["accepted_aphia_id"],
                    "accepted_name": resolved["accepted_name"],
                    "authority": resolved["authority"],
                    "rank": resolved["rank"],
                }
            )
            accepted_id = resolved["accepted_aphia_id"]
            bucket = taxa_rollup.setdefault(
                accepted_id,
                {
                    "accepted_aphia_id": accepted_id,
                    "accepted_name": resolved["accepted_name"],
                    "authority": resolved["authority"],
                    "rank": resolved["rank"],
                    "mention_count": 0,
                    "first_chunk": ch.get("chunk_id"),
                },
            )
            bucket["mention_count"] += 1

    taxa_list = sorted(
        taxa_rollup.values(),
        key=lambda r: (-r["mention_count"], r["accepted_name"]),
    )
    return {
        "total_mentions": len(mentions),
        "unique_taxa": len(taxa_list),
        "mentions": mentions,
        "taxa": taxa_list,
    }


# ---------------------------------------------------------------------------
# Anatomy lexicon and term extraction
# ---------------------------------------------------------------------------


def load_anatomy_lexicon(path: Path) -> Dict[str, Dict]:
    """Load the YAML lexicon into an in-memory dict.

    Structure returned:
        {canonical_term: {"synonyms": [...], "translations": {lang: [...]},
                          "description": str}}

    Keys are as-written in the YAML (lowercase snake_case recommended).
    """
    import yaml
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    out: Dict[str, Dict] = {}
    for canonical, entry in data.items():
        entry = entry or {}
        out[canonical] = {
            "synonyms": list(entry.get("synonyms") or []),
            "translations": dict(entry.get("translations") or {}),
            "description": entry.get("description") or "",
        }
    return out


def _build_anatomy_matcher(lexicon: Dict[str, Dict]) -> Tuple[re.Pattern, Dict[str, str]]:
    """Compile one big alternation regex and a variant→canonical map.

    Variants are matched as whole words, case-insensitive. Canonical
    names (e.g., ``nectophore``) ARE matched — their own keys count as
    variants of themselves so you don't have to repeat them in the
    ``synonyms`` list.
    """
    variant_to_canonical: Dict[str, str] = {}
    for canonical, entry in lexicon.items():
        # The canonical key itself (with underscores replaced by spaces
        # since the YAML convention is snake_case).
        canonical_display = canonical.replace("_", " ")
        variants = {canonical, canonical_display, *entry["synonyms"]}
        for lang_variants in entry.get("translations", {}).values():
            variants.update(lang_variants)
        for v in variants:
            if v:
                variant_to_canonical[v.lower()] = canonical

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


def extract_anatomy_mentions(
    chunks: List[Dict],
    lexicon: Dict[str, Dict],
) -> Dict:
    """Scan chunks for anatomy lexicon terms.

    Same output shape as :func:`extract_taxon_mentions` — flat mentions
    list plus a per-canonical-term rollup. Matches are case-insensitive
    and whole-word; see :func:`_build_anatomy_matcher`.
    """
    if not lexicon:
        return {"total_mentions": 0, "unique_terms": 0, "mentions": [], "terms": []}

    pattern, variant_to_canonical = _build_anatomy_matcher(lexicon)

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
    return {
        "total_mentions": len(mentions),
        "unique_terms": len(terms_list),
        "mentions": mentions,
        "terms": terms_list,
    }

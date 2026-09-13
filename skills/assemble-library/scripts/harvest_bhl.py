#!/usr/bin/env python3
"""Harvest Biodiversity Heritage Library **Part** records — the pre-1930
literature — into a local cache.

TEMPLATE — copied into a library repo by the `assemble-library` skill.

BHL is the only practical route to article-level records for older material, so
for a clade with a deep historical literature this script is the difference
between a corpus that starts in 1990 and one that starts in 1758.

Search terms live in `bhl_terms.yaml`, not here. They are clade-specific, and
they are **not the same as the OpenAlex queries**: BHL indexes what the older
literature was published under, so the segregate and historical genera matter
more here than the current name does. That is the whole reason this source
exists.

## Parts, not Items

A **Part** is an article-level segment — its own title, authors, page range,
often a minted DOI. An **Item** is a whole scanned volume. A 600-page flora that
mentions your genus once is an Item, and admitting it would put a book in the
library disguised as a paper. So Items are dropped by default; `--include-items`
keeps them when you actually want volumes.

## Catalogue vs full-text search, and why the answer depends on your clade

`searchtype=C` searches catalogue metadata; `searchtype=F` searches scanned
full text. Measured on the live API:

* For a **common** name, full-text is nearly all cost. The viburnum library
  measured 7,577 full-text hits across fourteen terms producing exactly *one*
  Part the catalogue search had not already found, because full-text matches
  every volume that mentions the word in passing.
* For a **rare** name it is worth having. `Bargmannia` returns 1 Part by
  catalogue and 9 hits by full-text, of which 4 are Parts.

So catalogue is the default and `--fulltext` is opt-in and page-capped rather
than exhausted. Decide by how common your term is, not by a blanket rule — and
prefer running full-text only for the rare historical names, which is where it
pays.

## Getting the scans is a separate, unsolved problem

This harvests **records**, not PDFs. BHL's `partpdf` and `itempdf` endpoints
return 403 to any non-interactive client, so these entries land in the bib with
a `bhl_part` id, a rights statement and no `file` field — want-list entries, by
design. Do not try to work around that; see `references/retrieval-ethics.md`.

Usage:

    source ~/.config/corpus/secrets.env      # or however BHL_API_KEY is set
    python harvest_bhl.py
    python harvest_bhl.py --fulltext         # add full-text, capped
    python harvest_bhl.py --dry-run          # count hits, fetch no metadata
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("pyyaml is required: it is in the library environment.yaml")

API = "https://www.biodiversitylibrary.org/api3"
REPO = Path(__file__).resolve().parents[1]
BUILD = REPO / "build"
CACHE = BUILD / "bhl_parts.jsonl"
STATE = BUILD / "bhl_state.json"

GAP = 0.35              # BHL asks for courtesy; this keeps well inside it
PAGE_SIZE = 100
CATALOGUE_PAGE_CAP = 200
FULLTEXT_PAGE_CAP = 10  # full-text does not terminate on a common term


def api(op: str, key: str, **params) -> list:
    """One BHL call. Returns the Result list, or [] on any failure."""
    params.update(op=op, apikey=key, format="json")
    url = f"{API}?{urllib.parse.urlencode(params)}"
    for attempt in range(4):
        time.sleep(GAP)
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "corpus-library/1.0"}
            )
            with urllib.request.urlopen(req, timeout=90) as resp:
                data = json.load(resp)
        except (urllib.error.URLError, json.JSONDecodeError) as exc:
            if attempt == 3:
                print(f"  ! {op} failed: {exc}", file=sys.stderr)
                return []
            time.sleep(3 * (attempt + 1))
            continue
        if data.get("Status") != "ok":
            # An invalid key surfaces here, not as an HTTP error. Say so
            # plainly: a silent empty harvest looks like "no such literature".
            print(f"  ! {op}: {data.get('ErrorMessage')}", file=sys.stderr)
            return []
        return data.get("Result") or []
    return []


def load_terms(path: Path) -> list[str]:
    if not path.exists():
        sys.exit(
            f"no {path.name} found.\n"
            "BHL search terms are clade-specific and belong in a config. They are\n"
            "NOT the OpenAlex queries — favour the segregate and historical genera\n"
            "the older literature published under. For example:\n\n"
            "  - Viburnum\n"
            "  - Oreinotinus\n"
            "  - Solenotinus\n"
        )
    spec = yaml.safe_load(path.read_text(encoding="utf-8")) or []
    terms = [str(t).strip() for t in spec if str(t).strip()]
    if not terms:
        sys.exit(f"{path.name} lists no terms")
    return terms


def search(term: str, searchtype: str, key: str) -> list[dict]:
    """Page one PublicationSearch to exhaustion, or to the cap."""
    cap = FULLTEXT_PAGE_CAP if searchtype == "F" else CATALOGUE_PAGE_CAP
    out: list[dict] = []
    page = 1
    while True:
        got = api("PublicationSearch", key, searchterm=term,
                  searchtype=searchtype, page=page, pageSize=PAGE_SIZE)
        if not got:
            break
        out += got
        if len(got) < PAGE_SIZE:
            break
        page += 1
        if page > cap:
            print(f"    (capped {term}/{searchtype} at {cap * PAGE_SIZE} hits)",
                  file=sys.stderr)
            break
    return out


def load_state() -> dict:
    return json.loads(STATE.read_text(encoding="utf-8")) if STATE.exists() else {}


def save_state(state: dict) -> None:
    BUILD.mkdir(exist_ok=True)
    STATE.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")


def cached_part_ids() -> set[str]:
    if not CACHE.exists():
        return set()
    ids = set()
    with CACHE.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    ids.add(str(json.loads(line)["PartID"]))
                except (json.JSONDecodeError, KeyError):
                    continue
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--terms", default=None, help="terms YAML (default: bhl_terms.yaml beside this script)")
    ap.add_argument("--fulltext", action="store_true", help="also run full-text search (capped)")
    ap.add_argument("--include-items", action="store_true", help="keep whole volumes, not just Parts")
    ap.add_argument("--dry-run", action="store_true", help="count hits; fetch no Part metadata")
    args = ap.parse_args()

    key = (os.environ.get("BHL_API_KEY") or "").strip()
    if not key:
        # Degrade, do not abort: the rest of the harvest is still worth running,
        # and an unreported skip is indistinguishable from an absent literature.
        print(
            "BHL_API_KEY is not set — SKIPPING the BHL harvest.\n"
            "  This is where the pre-1930 literature comes from; without it the\n"
            "  library will start at whatever the modern indexes cover.\n"
            "  Free key: https://www.biodiversitylibrary.org/getapikey.aspx\n"
            "  See references/retrieval-ethics.md for how to set it safely.",
            file=sys.stderr,
        )
        return 0

    tpath = Path(args.terms) if args.terms else Path(__file__).resolve().parent / "bhl_terms.yaml"
    terms = load_terms(tpath)
    searchtypes = ["C", "F"] if args.fulltext else ["C"]

    BUILD.mkdir(exist_ok=True)
    state = load_state()
    known = cached_part_ids()

    # --- Phase 1: search, collecting Part ids -------------------------------
    found: dict[str, str] = {}   # PartID -> the term that found it
    n_items_dropped = 0
    for term in terms:
        for st in searchtypes:
            hits = search(term, st, key)
            parts = 0
            for h in hits:
                if h.get("BHLType") != "Part" and not args.include_items:
                    n_items_dropped += 1
                    continue
                pid = str(h.get("PartID") or "").strip()
                if pid:
                    found.setdefault(pid, term)
                    parts += 1
            print(f"  {term}/{st}: {len(hits)} hits, {parts} parts", file=sys.stderr)

    new = [p for p in found if p not in known]
    print(
        f"\n{len(found):,} parts matched, {len(new):,} new "
        f"({len(known):,} already cached, {n_items_dropped:,} items dropped)",
        file=sys.stderr,
    )

    if args.dry_run:
        print("--dry-run: no Part metadata fetched", file=sys.stderr)
        return 0

    # --- Phase 2: full metadata for each new Part ---------------------------
    # One call per Part, so this is the expensive half — hence resumability.
    added = 0
    with CACHE.open("a", encoding="utf-8") as sink:
        for i, pid in enumerate(new, 1):
            got = api("GetPartMetadata", key, id=pid, pages="f", names="f")
            if not got:
                continue
            record = got[0] if isinstance(got, list) else got
            if not (record.get("Title") or "").strip():
                continue
            record["_term"] = found[pid]
            sink.write(json.dumps(record, ensure_ascii=False) + "\n")
            added += 1
            if i % 25 == 0:
                sink.flush()
                state["last_part"] = pid
                save_state(state)
                print(f"    {i}/{len(new)} fetched", file=sys.stderr)

    state["last_part"] = new[-1] if new else state.get("last_part")
    save_state(state)
    print(f"\n{added:,} part records written to {CACHE.relative_to(REPO)}", file=sys.stderr)
    print(
        "These are records, not PDFs: BHL refuses scan downloads to "
        "non-interactive clients, so they become want-list entries.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

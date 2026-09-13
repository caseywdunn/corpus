#!/usr/bin/env python3
"""Find open-access PDF URLs for bib entries that do not have one yet, and write
the fetch plan `fetch_pdfs.py` reads.

TEMPLATE — copied into a library repo by the `assemble-library` skill.

OpenAlex's `best_oa_location` covers roughly a quarter of a typical bib.
Unpaywall indexes a wider, overlapping slice — in particular the "bronze" copies
publishers post without a licence, which OpenAlex's snapshot often lacks. Both
are keyless; both want `CORPUS_CONTACT_EMAIL` and give better limits for it.

This is deliberately separate from `fetch_pdfs.py`. Resolving is cheap, idempotent
and safe to re-run as the bib grows; fetching is slow, rude to repeat, and the
thing you want to be able to interrupt. Keeping them apart means a resolver
improving does not force a re-download, and a failed fetch does not lose the
resolution work.

Nothing here downloads a PDF. It writes candidate URLs and lets
`fetch_pdfs.py` decide whether what comes back is really a paper.

Usage:

    export CORPUS_CONTACT_EMAIL="you@example.edu"
    python resolve_oa.py                      # resolve everything unresolved
    python resolve_oa.py --limit 20           # a sample first
    python resolve_oa.py --refresh            # re-resolve, e.g. after a year
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

REPO = Path(__file__).resolve().parents[1]
BUILD = REPO / "build"
PLAN = BUILD / "fetch_plan.json"
CACHE = BUILD / "oa_resolutions.json"

UNPAYWALL = "https://api.unpaywall.org/v2/"
GAP = 0.2


def contact() -> str:
    value = (os.environ.get("CORPUS_CONTACT_EMAIL") or "").strip()
    if not value:
        sys.exit(
            "CORPUS_CONTACT_EMAIL is not set. Unpaywall requires it, and every "
            "other service here rewards it.\n"
            '  export CORPUS_CONTACT_EMAIL="you@example.edu"'
        )
    return value


def unpaywall_pdf_urls(doi: str, email: str) -> list[str]:
    """OA PDF locations Unpaywall knows for one DOI, best first."""
    url = f"{UNPAYWALL}{urllib.parse.quote(doi)}?{urllib.parse.urlencode({'email': email})}"
    req = urllib.request.Request(url, headers={"User-Agent": f"corpus-library/1.0 (mailto:{email})"})
    try:
        with urllib.request.urlopen(req, timeout=45) as resp:
            data = json.load(resp)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return []          # not in Unpaywall; ordinary, not an error
        if exc.code == 422:
            return []          # malformed DOI; the bib's problem, not ours
        raise
    except (urllib.error.URLError, json.JSONDecodeError):
        return []

    out: list[str] = []
    best = data.get("best_oa_location") or {}
    if best.get("url_for_pdf"):
        out.append(best["url_for_pdf"])
    for loc in data.get("oa_locations") or []:
        pdf = loc.get("url_for_pdf")
        if pdf and pdf not in out:
            out.append(pdf)
    return out


def load_json(path: Path, default):
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else default


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--entries", default=None,
                    help="JSON {filename: {doi, oa_url?}} to resolve; default build/bib_entries.json")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--refresh", action="store_true", help="ignore cached resolutions")
    args = ap.parse_args()

    email = contact()
    entries_path = Path(args.entries) if args.entries else BUILD / "bib_entries.json"
    if not entries_path.exists():
        sys.exit(
            f"no {entries_path.name}. Emit it from the bib first — the entries\n"
            "that have no PDF yet, each with whatever identifiers are known:\n\n"
            '  {"Smith1998.pdf": {"doi": "10.1234/x", "oa_url": null}}\n'
        )
    entries = json.loads(entries_path.read_text(encoding="utf-8"))
    cache = {} if args.refresh else load_json(CACHE, {})
    plan = load_json(PLAN, {})

    resolved = already = nothing = 0
    for i, (fname, row) in enumerate(sorted(entries.items())):
        if args.limit and i >= args.limit:
            break

        # An OA URL the harvest already knew beats a lookup.
        candidates: list[str] = []
        if row.get("oa_url"):
            candidates.append(row["oa_url"])

        doi = (row.get("doi") or "").strip()
        if doi:
            if doi in cache and not args.refresh:
                candidates += [u for u in cache[doi] if u not in candidates]
                already += 1
            else:
                time.sleep(GAP)
                try:
                    found = unpaywall_pdf_urls(doi, email)
                except urllib.error.HTTPError as exc:
                    print(f"  ! {doi}: HTTP {exc.code}", file=sys.stderr)
                    found = []
                cache[doi] = found
                candidates += [u for u in found if u not in candidates]
                if found:
                    resolved += 1

        if candidates:
            # fetch_pdfs.py takes one URL; keep the rest so a later pass can
            # try the next one without re-resolving.
            plan[fname] = {"url": candidates[0], "doi": doi or None,
                           "alternates": candidates[1:]}
        else:
            nothing += 1

    BUILD.mkdir(exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")
    PLAN.write_text(json.dumps(plan, indent=2, sort_keys=True), encoding="utf-8")

    print(
        f"{len(plan):,} entries in the fetch plan "
        f"({resolved:,} newly resolved, {already:,} cached, {nothing:,} with no OA copy)",
        file=sys.stderr,
    )
    print(
        "Entries with no OA copy are want-list entries, not failures — most "
        "literature is not open access, and the bib keeps them either way.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

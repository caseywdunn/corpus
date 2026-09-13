#!/usr/bin/env python3
"""Harvest OpenAlex into a local record cache — incremental, resumable, and
aware of what it is spending.

TEMPLATE — copied into a library repo by the `assemble-library` skill.

The queries are **not** in this file. They are clade-specific and they are the
output of the scoping conversation, so they live in `queries.yaml` beside it
(see `--queries`). That separation is the point: someone revising the harvest
for a new clade edits a list of queries, not a program.

## What this knows about OpenAlex that a naive harvest does not

OpenAlex meters by **credit**, not by call, and the free allowance is small —
about 1000 credits / $0.10 a day, where a filtered search costs 10 credits and a
plain listing costs 1. That is roughly 100 search requests per day. Run
`openalex_probe.py --check` for the live numbers and `--estimate` to cost a
query before adding it here.

Three consequences are built in:

* **`per_page=200` always.** Cost is per request, not per record, so a smaller
  page multiplies the bill for identical results.
* **A 429 is not always transient.** When the daily allowance is spent, the
  `Retry-After` is measured in *hours*. Retrying through that is pointless and
  rude, so this stops cleanly and leaves resumable state instead. A short
  `Retry-After` — the per-second limiter — is waited out normally.
* **Resume is per query and per page.** State records the cursor each query
  reached, so tomorrow's run continues from there rather than re-paging ground
  it already paid for.

`select=` keeps the payload to the fields actually parsed. It does not change
the price, but it does change how long a large harvest takes.

Usage:

    export CORPUS_CONTACT_EMAIL="you@example.edu"

    python harvest_openalex.py                    # harvest / resume
    python harvest_openalex.py --dry-run          # cost it first, spend nothing
    python harvest_openalex.py --max-requests 40  # stay inside a budget
    python harvest_openalex.py --refresh          # discard state, start over
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

API = "https://api.openalex.org/works"
REPO = Path(__file__).resolve().parents[1]
BUILD = REPO / "build"
CACHE = BUILD / "records.jsonl"
STATE = BUILD / "harvest_state.json"

PER_PAGE = 200          # the maximum, and the only sensible value — see above
PAGE_DELAY = 0.5        # the polite pool allows ~10 req/s; sustained paging needs less
SHORT_RETRY_CEILING = 120   # a Retry-After above this means "come back tomorrow"

SELECT = ",".join([
    "id", "doi", "ids", "display_name", "title", "publication_year",
    "publication_date", "type", "language", "authorships", "biblio",
    "primary_location", "best_oa_location", "open_access", "is_retracted",
    "referenced_works",
])


def contact() -> str:
    value = (os.environ.get("CORPUS_CONTACT_EMAIL") or "").strip()
    if not value:
        sys.exit(
            "CORPUS_CONTACT_EMAIL is not set. Identified clients get materially "
            "higher rate limits.\n"
            '  export CORPUS_CONTACT_EMAIL="you@example.edu"'
        )
    return value


class Exhausted(Exception):
    """The daily allowance is gone. Not an error — a place to stop."""

    def __init__(self, retry_after: int | None):
        self.retry_after = retry_after
        super().__init__("daily allowance exhausted")


def get(params: dict, email: str) -> tuple[dict, dict]:
    """One request, with the retry policy the credit model actually calls for."""
    query = dict(params)
    query["mailto"] = email
    url = f"{API}?{urllib.parse.urlencode(query, safe=':,%\"')}"
    req = urllib.request.Request(
        url, headers={"User-Agent": f"corpus-library/1.0 (mailto:{email})"}
    )
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                headers = {k.lower(): v for k, v in resp.headers.items()}
                return json.load(resp), headers
        except urllib.error.HTTPError as exc:
            if exc.code == 429:
                raw = exc.headers.get("Retry-After")
                try:
                    wait = int(raw)
                except (TypeError, ValueError):
                    wait = None
                # The distinction that matters: a short wait is the per-second
                # limiter and is worth sleeping through; a long one means the
                # daily credit allowance is spent, and no amount of retrying
                # will change that before it resets.
                if wait is not None and wait <= SHORT_RETRY_CEILING:
                    time.sleep(wait + 1)
                    continue
                raise Exhausted(wait) from exc
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))
        except urllib.error.URLError:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("unreachable")


def budget(headers: dict) -> tuple[float | None, float | None]:
    def num(key):
        try:
            return float(headers.get(key))
        except (TypeError, ValueError):
            return None
    return num("x-ratelimit-remaining-usd"), num("x-ratelimit-cost-usd")


def load_queries(path: Path) -> list[tuple[str, str]]:
    """`queries.yaml` — a list of {tag, filter}, written per clade.

    The tag travels with every record it found, because the relevance filter
    downstream should be able to reason over *how* something was found: a
    title-search hit is a different kind of evidence from a vernacular sweep,
    and a sweep that turns out to be contaminated can be quarantined by tag
    rather than by rebuilding the cache.
    """
    if not path.exists():
        sys.exit(
            f"no {path.name} found.\n"
            "The queries are clade-specific and belong in a config, not in this "
            "script. Write one from the scoping conversation, e.g.:\n\n"
            "  - tag: genus-title\n"
            "    filter: title.search:Viburnum\n"
            "  - tag: vernacular-guelder\n"
            '    filter: title.search:%22guelder rose%22\n'
        )
    spec = yaml.safe_load(path.read_text(encoding="utf-8")) or []
    out = []
    for item in spec:
        if not isinstance(item, dict) or "tag" not in item or "filter" not in item:
            sys.exit(f"malformed entry in {path.name}: {item!r} (need tag + filter)")
        out.append((str(item["tag"]), str(item["filter"])))
    if not out:
        sys.exit(f"{path.name} lists no queries")
    return out


def load_state() -> dict:
    if STATE.exists():
        return json.loads(STATE.read_text(encoding="utf-8"))
    return {}


def save_state(state: dict) -> None:
    BUILD.mkdir(exist_ok=True)
    STATE.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")


def seen_ids() -> set[str]:
    """Work ids already cached, so a resumed run does not duplicate them."""
    if not CACHE.exists():
        return set()
    ids = set()
    with CACHE.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    ids.add(json.loads(line)["id"])
                except (json.JSONDecodeError, KeyError):
                    continue
    return ids


def dry_run(queries, email) -> int:
    """Cost the whole harvest without paging any of it. One request per query."""
    total_pages = 0
    total_cost = 0.0
    print(f"{'tag':<32} {'matches':>10} {'requests':>9}")
    for tag, filt in queries:
        payload, headers = get({"filter": filt, "per-page": 1}, email)
        count = (payload.get("meta") or {}).get("count", 0)
        pages = -(-count // PER_PAGE)
        cost = (payload.get("meta") or {}).get("cost_usd") or 0.0
        total_pages += pages
        total_cost += pages * cost
        print(f"{tag:<32} {count:>10,} {pages:>9,}")
        time.sleep(PAGE_DELAY)
    remaining, _ = budget(headers)
    print(f"\n{'TOTAL':<32} {'':>10} {total_pages:>9,} requests")
    if total_cost:
        print(f"estimated cost: ${total_cost:.4f}")
        if remaining is not None:
            print(f"remaining today: ${remaining:.4f}")
            if total_cost > remaining:
                days = total_cost / 0.10
                print(
                    f"This harvest needs ~{days:.1f} days of free allowance. "
                    "It will resume across runs; that is expected, not a failure."
                )
    return 0


def harvest(queries, email, max_requests: int | None) -> int:
    BUILD.mkdir(exist_ok=True)
    state = load_state()
    known = seen_ids()
    requests_made = 0
    added = 0
    stopped_early = False

    with CACHE.open("a", encoding="utf-8") as sink:
        for tag, filt in queries:
            q = state.setdefault(tag, {"done": False, "cursor": "*", "seen": 0})
            if q["done"]:
                print(f"  {tag}: done", file=sys.stderr)
                continue

            while q["cursor"]:
                if max_requests is not None and requests_made >= max_requests:
                    print(f"\nreached --max-requests {max_requests}", file=sys.stderr)
                    stopped_early = True
                    break
                time.sleep(PAGE_DELAY)
                try:
                    page, headers = get({
                        "filter": filt,
                        "select": SELECT,
                        "per-page": PER_PAGE,
                        "cursor": q["cursor"],
                    }, email)
                except Exhausted as exc:
                    hours = f"{exc.retry_after / 3600:.1f} h" if exc.retry_after else "unknown"
                    print(
                        f"\nDaily OpenAlex allowance spent (Retry-After: {hours}).\n"
                        "State saved — re-run tomorrow and the harvest continues "
                        "from here rather than restarting.",
                        file=sys.stderr,
                    )
                    stopped_early = True
                    break
                requests_made += 1

                results = page.get("results") or []
                for record in results:
                    if record.get("id") in known:
                        continue
                    known.add(record["id"])
                    record["_query"] = tag
                    sink.write(json.dumps(record, ensure_ascii=False) + "\n")
                    added += 1

                q["seen"] += len(results)
                q["cursor"] = (page.get("meta") or {}).get("next_cursor")
                if not results or not q["cursor"]:
                    q["done"] = True
                    q["cursor"] = None
                total = (page.get("meta") or {}).get("count", 0)
                print(f"    {tag}: {q['seen']:,}/{total:,}", file=sys.stderr)
                sink.flush()
                save_state(state)

            save_state(state)
            if stopped_early:
                break

    save_state(state)
    done = sum(1 for v in state.values() if v.get("done"))
    print(
        f"\n{added:,} new records ({len(known):,} cached), "
        f"{requests_made} requests, {done}/{len(queries)} queries complete",
        file=sys.stderr,
    )
    if stopped_early:
        print("Incomplete — re-run to continue.", file=sys.stderr)
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--queries", default=None, help="queries YAML (default: queries.yaml beside this script)")
    ap.add_argument("--dry-run", action="store_true", help="cost the harvest without paging it")
    ap.add_argument("--max-requests", type=int, default=None, help="stop after this many requests")
    ap.add_argument("--refresh", action="store_true", help="discard resume state and start over")
    args = ap.parse_args()

    email = contact()
    qpath = Path(args.queries) if args.queries else Path(__file__).resolve().parent / "queries.yaml"
    queries = load_queries(qpath)

    if args.refresh and STATE.exists():
        STATE.unlink()
        print("resume state discarded (records.jsonl kept; it deduplicates)", file=sys.stderr)

    if args.dry_run:
        return dry_run(queries, email)
    return harvest(queries, email, args.max_requests)


if __name__ == "__main__":
    raise SystemExit(main())

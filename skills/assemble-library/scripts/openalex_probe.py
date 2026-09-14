#!/usr/bin/env python3
"""Probe OpenAlex cheaply: is it reachable, what does a harvest cost, does the
response still carry the fields a harvest depends on?

TEMPLATE — copied into a library repo by the `assemble-library` skill.

Run this **before** writing or running a harvest. It exists because OpenAlex
meters by credit rather than by call, the free allowance is small, and the
difference between a well-planned harvest and a careless one is the difference
between finishing today and waiting until tomorrow.

Measured against the live API (2026-09-13), with the numbers the headers
themselves report:

    x-ratelimit-limit          1000 credits/day
    x-ratelimit-limit-usd      0.10 USD/day
    plain listing              0.0001 USD   (1 credit)
    filtered / search query    0.001  USD   (10 credits)

So the free tier is roughly **100 search requests per day**, not the hundreds of
thousands of plain calls an older reading of the docs would suggest. Two
consequences the harvest must respect:

* **Page size is the whole game.** `per_page=200` is the maximum and costs the
  same as `per_page=25`, so a small page size throws away 8x your daily budget
  for nothing.
* **Never page blindly.** `meta.count` tells you how many works a query matches
  from a *single* request, so the cost of a harvest is knowable before you spend
  it. `--estimate` does exactly that, and it is the reason this script exists.

Reported costs are what the API returned at the time of measurement; treat them
as the current shape of the model rather than a permanent constant, which is why
this script reads the live headers instead of hardcoding them.

Usage:

    export CORPUS_CONTACT_EMAIL="you@example.edu"

    python openalex_probe.py --check
    python openalex_probe.py --estimate 'title.search:Viburnum'
    python openalex_probe.py --estimate 'title.search:Viburnum' --per-page 200
    python openalex_probe.py --contract          # do we still get the fields we parse?

Exits non-zero on a failed check so it can gate a harvest in a script.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.openalex.org/works"

# Fields the harvest and the reference-mining step actually read. Probing these
# is cheaper than discovering mid-harvest that a rename broke the parse.
REQUIRED_FIELDS = (
    "id",
    "doi",
    "display_name",
    "publication_year",
    "type",
    "authorships",
    "referenced_works",
    # A harvest that omits this silently defeats any abstract-based relevance
    # rule, so the contract check asserts it is still served.
    "abstract_inverted_index",
)

# The maximum OpenAlex accepts, and the only page size worth using: cost is per
# request, not per record.
MAX_PER_PAGE = 200


def contact() -> str:
    """The polite-pool address. Failing here is deliberate — see below."""
    value = (os.environ.get("CORPUS_CONTACT_EMAIL") or "").strip()
    if not value:
        sys.exit(
            "CORPUS_CONTACT_EMAIL is not set. Public indexes give identified "
            "clients materially higher rate limits, and an unidentified harvest "
            "is both slower and ruder.\n"
            '  export CORPUS_CONTACT_EMAIL="you@example.edu"\n'
            "There is no safe default: a placeholder lies to the service, and "
            "someone else's address sends them the consequences of your traffic."
        )
    return value


def request(params: dict) -> tuple[dict, dict]:
    """One GET. Returns (payload, headers). Never retries — a probe that
    hammers a service it is checking the politeness of has missed the point."""
    email = contact()
    query = dict(params)
    query["mailto"] = email
    url = f"{API}?{urllib.parse.urlencode(query)}"
    req = urllib.request.Request(
        url, headers={"User-Agent": f"corpus-library-probe/1.0 (mailto:{email})"}
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            # Lowercase the keys: HTTP header names are case-insensitive, but a
            # plain dict lookup is not, and OpenAlex sends them in mixed case.
            headers = {k.lower(): v for k, v in resp.headers.items()}
            return json.load(resp), headers
    except urllib.error.HTTPError as exc:
        if exc.code == 429:
            retry = exc.headers.get("Retry-After", "unknown")  # Message: case-insensitive
            sys.exit(
                f"429 from OpenAlex; Retry-After: {retry}s. The daily credit "
                "allowance is spent. A harvest must be incremental and resumable "
                "so tomorrow's run continues rather than restarts."
            )
        sys.exit(f"HTTP {exc.code} from OpenAlex: {exc.reason}")
    except urllib.error.URLError as exc:
        sys.exit(f"could not reach OpenAlex: {exc.reason}")


def _budget(headers: dict) -> dict:
    """Pull the credit picture out of the response headers."""
    def num(key, cast=float):
        raw = headers.get(key)
        try:
            return cast(raw)
        except (TypeError, ValueError):
            return None

    return {
        "credits_limit": num("x-ratelimit-limit", int),
        "credits_remaining": num("x-ratelimit-remaining", int),
        "usd_limit": num("x-ratelimit-limit-usd"),
        "usd_remaining": num("x-ratelimit-remaining-usd"),
        "this_request_usd": num("x-ratelimit-cost-usd"),
        "reset_seconds": num("x-ratelimit-reset", int),
    }


def _report_budget(b: dict) -> None:
    if b["credits_limit"]:
        print(
            f"  credits:   {b['credits_remaining']}/{b['credits_limit']} remaining"
        )
    if b["usd_limit"]:
        print(
            f"  allowance: ${b['usd_remaining']:.4f} of ${b['usd_limit']:.2f} remaining"
        )
    if b["this_request_usd"] is not None:
        print(f"  this call: ${b['this_request_usd']:.4f}")
    if b["reset_seconds"]:
        hours = b["reset_seconds"] / 3600
        print(f"  resets in: {hours:.1f} h")


def cmd_check() -> int:
    """Reachable, identified, and how much budget is left. One request."""
    payload, headers = request({"per-page": 1})
    print("OpenAlex reachable, polite pool accepted.")
    budget = _budget(headers)
    _report_budget(budget)

    if budget["usd_remaining"] is not None and budget["this_request_usd"]:
        # A search costs ~10x a plain listing, so quote the honest figure.
        search_cost = budget["this_request_usd"] * 10
        affordable = int(budget["usd_remaining"] // search_cost)
        print(
            f"\n  ≈{affordable} search requests left today "
            f"(a filtered query costs ~${search_cost:.4f}, "
            f"{search_cost / budget['this_request_usd']:.0f}x a plain listing)"
        )
        if affordable < 20:
            print("  Low. Plan the harvest, or resume tomorrow.")
    return 0


def cmd_estimate(query: str, per_page: int) -> int:
    """What a harvest of this query would cost — from ONE request.

    This is the function that justifies the script: `meta.count` makes the size
    of a harvest knowable before committing to it.
    """
    payload, headers = request({"filter": query, "per-page": 1})
    meta = payload.get("meta") or {}
    count = meta.get("count")
    if count is None:
        sys.exit(f"no meta.count in response for {query!r}; check the filter syntax")

    pages = math.ceil(count / per_page) if count else 0
    per_call = meta.get("cost_usd")
    budget = _budget(headers)
    if per_call is None:
        per_call = budget["this_request_usd"]

    print(f"query:   {query}")
    print(f"matches: {count:,} works")
    print(f"harvest: {pages:,} requests at per_page={per_page}")
    if per_call:
        total = pages * per_call
        print(f"cost:    ~${total:.4f} (${per_call:.4f} per request)")
        if budget["usd_limit"]:
            days = total / budget["usd_limit"]
            if days > 1:
                print(
                    f"         ~{days:.1f} days of the free allowance — "
                    "the harvest must be resumable"
                )
    if per_page < MAX_PER_PAGE:
        wasted = math.ceil(count / per_page) - math.ceil(count / MAX_PER_PAGE)
        print(
            f"\n  per_page={per_page} costs {wasted:,} more requests than "
            f"per_page={MAX_PER_PAGE}, for identical results."
        )
    _report_budget(budget)
    return 0


def cmd_contract() -> int:
    """Does a work record still carry the fields the harvest parses?

    Cheap insurance: a renamed field surfaces here as one failed check rather
    than as a harvest that runs to completion and writes empty columns.
    """
    payload, headers = request({"per-page": 1})
    results = payload.get("results") or []
    if not results:
        sys.exit("no results returned; cannot check the field contract")
    work = results[0]

    missing = [f for f in REQUIRED_FIELDS if f not in work]
    for field in REQUIRED_FIELDS:
        print(f"  {'ok  ' if field in work else 'MISS'}  {field}")

    # Nested shapes the harvest reaches into, checked only where present so a
    # record that legitimately lacks one is not reported as a broken contract.
    authorships = work.get("authorships") or []
    if authorships and "author" not in authorships[0]:
        missing.append("authorships[].author")
        print("  MISS  authorships[].author")
    elif authorships:
        print("  ok    authorships[].author")

    _report_budget(_budget(headers))
    if missing:
        print(f"\nFAILED: {len(missing)} field(s) missing: {', '.join(missing)}")
        print("The harvest parses these; fix the parse before running it.")
        return 1
    print("\nField contract holds.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="reachability + remaining budget")
    ap.add_argument("--estimate", metavar="FILTER", help="cost a harvest of this filter")
    ap.add_argument("--contract", action="store_true", help="verify parsed fields still exist")
    ap.add_argument(
        "--per-page", type=int, default=MAX_PER_PAGE,
        help=f"page size for the estimate (default and maximum {MAX_PER_PAGE})",
    )
    args = ap.parse_args()

    if args.per_page > MAX_PER_PAGE:
        sys.exit(f"--per-page cannot exceed {MAX_PER_PAGE}")

    if not (args.check or args.estimate or args.contract):
        ap.print_help()
        return 2

    rc = 0
    if args.check:
        rc |= cmd_check()
    if args.estimate:
        rc |= cmd_estimate(args.estimate, args.per_page)
    if args.contract:
        rc |= cmd_contract()
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Fetch open-access PDFs for bib entries that do not have one yet.

TEMPLATE — copied into a library repo by the `assemble-library` skill.

Most attempts fail, and that is the designed outcome, not a broken run. Expect
roughly a quarter to a third to succeed. Every failure becomes a **want-list**
entry: the bib keeps the paper with its DOI or URL and no `file` field, and a
human fetches it later with institutional access. See
`references/retrieval-ethics.md` — in particular, a 403 is a final answer and
BHL scan endpoints are deliberately not attempted.

## Validation is the point of this script

A download that returns bytes is not a paper. Publisher paywalls return HTTP
200 with an HTML interstitial, and those must never land in `library/` looking
like literature — a corpus build would ingest one as the paper its bib entry
names, and nothing downstream would notice. So a candidate is kept only if it:

1. starts with `%PDF-`,
2. clears a floor size (a few kB of "access denied" is not an article), and
3. parses under `pdfinfo`.

**Know what this does and does not catch**, measured rather than assumed:

* Caught — HTML interstitials, stubs below the floor, non-PDF bytes, truncated
  downloads, a missing trailer or a broken xref. That covers the realistic
  network failures.
* **Not caught** — a structurally valid PDF of the *wrong paper*. `pdfinfo` is
  lenient about semantic damage: changing an object's `/Type` passes cleanly as
  long as the xref still resolves. So this guarantees "a readable PDF arrived",
  never "the right PDF arrived".

That second gap is real and has bitten a library before — a mis-resolved
reference lands a genuine PDF of a different paper under a plausible filename,
and nothing downstream notices. It is closed after a build by comparing each
bib title against the PDF's extracted text; see
`references/bib-conventions.md`, "Verifying a binding actually holds".

`is_acceptable_pdf()` is a pure function over bytes so it can be tested without
the network, which is most of what is worth testing here.

Usage:

    python fetch_pdfs.py --dry-run          # what would be attempted
    python fetch_pdfs.py --limit 5          # fetch a few
    python fetch_pdfs.py                    # everything outstanding
    python fetch_pdfs.py --retry-failed     # after a resolver improves
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LIBRARY = REPO / "library"
BUILD = REPO / "build"
FAILURES = BUILD / "fetch_failures.json"

MIN_BYTES = 20_000      # smaller than this is a stub, a cover page, or an error
TIMEOUT = 120
GAP = 1.0               # per-host courtesy pause

# Endpoints known to refuse non-interactive clients. Attempting them only wastes
# requests and looks like probing; see references/retrieval-ethics.md.
NEVER_ATTEMPT = ("biodiversitylibrary.org/partpdf", "biodiversitylibrary.org/itempdf")


def is_acceptable_pdf(blob: bytes, path: Path | None = None) -> tuple[bool, str]:
    """Is this actually a usable PDF? Returns (ok, reason).

    Pure over `blob` except for the optional `pdfinfo` structural check, which
    needs a file on disk. Kept separate from any network code precisely so the
    rejection paths can be exercised with fabricated inputs.
    """
    if len(blob) < MIN_BYTES:
        return False, f"too small ({len(blob):,} bytes < {MIN_BYTES:,})"

    head = blob[:1024].lstrip()
    if not (blob[:5] == b"%PDF-" or head[:5] == b"%PDF-"):
        # Name the common case rather than saying "not a PDF": an operator
        # seeing "html interstitial" knows it is a paywall, not a broken link.
        lowered = head[:200].lower()
        if b"<html" in lowered or b"<!doctype" in lowered:
            return False, "html interstitial (paywall or login page)"
        return False, "not a PDF (no %PDF- header)"

    if path is not None:
        try:
            probe = subprocess.run(
                ["pdfinfo", str(path)], capture_output=True, timeout=60
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            return False, f"pdfinfo failed to run: {exc}"
        if probe.returncode != 0:
            detail = probe.stderr.decode("utf-8", "replace").strip().splitlines()
            return False, f"pdfinfo rejected it: {detail[0] if detail else 'unknown'}"

    return True, "ok"


def download(url: str, dest: Path, contact: str) -> tuple[bool, str]:
    """Fetch one candidate URL into `dest`. Never retries a refusal."""
    if any(marker in url for marker in NEVER_ATTEMPT):
        return False, "endpoint refuses non-interactive clients (not attempted)"

    req = urllib.request.Request(
        url,
        headers={
            # A real, identifying User-Agent. Never spoofed to evade a block —
            # if a host refuses this, that is its answer.
            "User-Agent": f"corpus-library/1.0 (+mailto:{contact})",
            "Accept": "application/pdf,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            blob = resp.read()
    except urllib.error.HTTPError as exc:
        if exc.code in (401, 402, 403):
            return False, f"HTTP {exc.code} (access refused — final)"
        return False, f"HTTP {exc.code}"
    except (urllib.error.URLError, TimeoutError) as exc:
        return False, f"network error: {exc}"

    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        tmp.write(blob)
        tmp_path = Path(tmp.name)
    try:
        ok, reason = is_acceptable_pdf(blob, tmp_path)
        if ok:
            shutil.move(str(tmp_path), dest)
            return True, f"{len(blob):,} bytes"
        return False, reason
    finally:
        tmp_path.unlink(missing_ok=True)


def shelf_for(stem: str) -> Path:
    """`library/<letter>/` — keeps a listing usable at a few thousand files."""
    first = next((c for c in stem if c.isalpha()), None)
    return LIBRARY / (first.upper() if first else "_")


def load_failures() -> dict:
    return json.loads(FAILURES.read_text(encoding="utf-8")) if FAILURES.exists() else {}


def save_failures(failures: dict) -> None:
    BUILD.mkdir(exist_ok=True)
    FAILURES.write_text(json.dumps(failures, indent=2, sort_keys=True), encoding="utf-8")


def main() -> int:
    import os

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plan", default=None,
                    help="JSON mapping {filename: {url: ...}}; default build/fetch_plan.json")
    ap.add_argument("--limit", type=int, default=None, help="stop after N successful fetches")
    ap.add_argument("--dry-run", action="store_true", help="list what would be attempted")
    ap.add_argument("--retry-failed", action="store_true",
                    help="re-attempt cached failures (use after a resolver improves)")
    args = ap.parse_args()

    contact = (os.environ.get("CORPUS_CONTACT_EMAIL") or "").strip()
    if not contact:
        sys.exit(
            'CORPUS_CONTACT_EMAIL is not set — it identifies this client to hosts.\n'
            '  export CORPUS_CONTACT_EMAIL="you@example.edu"'
        )

    plan_path = Path(args.plan) if args.plan else BUILD / "fetch_plan.json"
    if not plan_path.exists():
        sys.exit(
            f"no {plan_path.name}. Build one from the bib first — it maps the\n"
            "target filename to the candidate URLs resolved for that entry:\n\n"
            '  {"Smith1998.pdf": {"url": "https://...", "doi": "10.1234/x"}}\n'
        )
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    failures = {} if args.retry_failed else load_failures()

    attempted = fetched = skipped = 0
    for fname, row in sorted(plan.items()):
        dest = shelf_for(fname) / fname
        if dest.exists() and dest.stat().st_size >= MIN_BYTES:
            skipped += 1
            continue
        if fname in failures and not args.retry_failed:
            skipped += 1
            continue
        url = (row.get("url") or "").strip()
        if not url:
            failures[fname] = "no candidate URL"
            continue

        if args.dry_run:
            print(f"  would fetch {fname} <- {url[:90]}")
            attempted += 1
            continue

        time.sleep(GAP)
        ok, detail = download(url, dest, contact)
        attempted += 1
        if ok:
            fetched += 1
            failures.pop(fname, None)
            print(f"  ok    {fname}  ({detail})")
        else:
            failures[fname] = detail
            print(f"  fail  {fname}  {detail}")
        if args.limit and fetched >= args.limit:
            print(f"  (stopping at --limit {args.limit})")
            break

    if not args.dry_run:
        save_failures(failures)
        rate = (fetched / attempted * 100) if attempted else 0.0
        print(
            f"\n{fetched}/{attempted} fetched ({rate:.0f}%), {skipped} skipped, "
            f"{len(failures)} on the want-list"
        )
        print(
            "A failure is a want-list entry, not a defeat: the bib keeps the "
            "paper with no `file` field, and a human can fetch it later."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

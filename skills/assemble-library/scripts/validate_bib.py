#!/usr/bin/env python3
"""Validate a library's ``.bib`` against ``library/`` and surface issues.

TEMPLATE — copied into a library repo by the `assemble-library` skill. It is
clade-agnostic by construction: the bib is discovered rather than named, so
nothing here needs editing for a new group. Edit it freely once copied; it is
yours, not the tool's.

This is the day-to-day curation tool now that the bib is the source of
truth. The docx / orphan-reconciliation pipeline is for *initial* import;
once you've decided what's in the bib, this script tells you whether the
bib is internally consistent and inventoried correctly against the PDFs.

Output: ``logs/validate_bib.log`` (also printed to stdout). Sections:

    SUMMARY            high-level coverage numbers
    INVENTORY GAPS     bib ↔ library/ mismatches that need fixing
    METADATA HEALTH    missing/conflicting fields, duplicate keys/DOIs
    DATA INTEGRITY     malformed DOIs, URLs, etc.
    STATS BY DECADE    decade-by-decade coverage table
    PDF CONTENT SCAN   (only with --scan-pdfs) duplicate hashes, page
                       stats, corrupt files

Goal: a clean run is one where INVENTORY GAPS is empty and you've
explicitly accepted whatever METADATA HEALTH issues remain.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import bibtexparser

REPO = Path(__file__).resolve().parents[1]
LIBRARY = REPO / "library"
LOGS = REPO / "logs"
LOGS.mkdir(exist_ok=True)
def _discover_bib(explicit: str | None = None) -> Path:
    """The library's .bib — named explicitly, else the only one at the root.

    Discovery rather than a constant because this file is a template: a
    hardcoded name is the first thing every copy would have to edit, and the
    edit that gets forgotten. Ambiguity is an error rather than a guess, since
    silently validating the wrong bib reports a clean run against the wrong
    inventory.
    """
    if explicit:
        path = Path(explicit)
        if not path.is_absolute():
            path = REPO / path
        if not path.exists():
            raise SystemExit(f"--bib {path} does not exist")
        return path
    candidates = sorted(p for p in REPO.glob("*.bib") if p.is_file())
    if not candidates:
        raise SystemExit(f"no .bib found at the root of {REPO}; pass --bib")
    if len(candidates) > 1:
        names = ", ".join(p.name for p in candidates)
        raise SystemExit(f"several .bib files at {REPO} ({names}); pass --bib")
    return candidates[0]


BIB = _discover_bib()  # rebound in main() once --bib is parsed
LIBRARY_NAME = REPO.name
OUT = LOGS / "validate_bib.log"
README = REPO / "readme.md"
ASSETS = REPO / "assets"
HISTOGRAM_PNG = ASSETS / "library_stats.png"
README_BEGIN = "<!-- BEGIN: stats (autogen by scripts/validate_bib.py --emit-readme) -->"
README_END = "<!-- END: stats -->"

DOI_RE = re.compile(r"^10\.\d{4,9}/.+$")
URL_RE = re.compile(r"^https?://[^\s]+$")


def collect_library_pdfs() -> set[str]:
    """All PDF basenames in library/, excluding the orphans/ subdir."""
    out: set[str] = set()
    for p in LIBRARY.rglob("*.pdf"):
        rel = p.relative_to(LIBRARY)
        if rel.parts and rel.parts[0] == "orphans":
            continue
        out.add(p.name)
    return out


def pdf_paths_by_basename() -> dict[str, list[Path]]:
    """basename → list of full paths (catches duplicate basenames in subdirs)."""
    out: dict[str, list[Path]] = defaultdict(list)
    for p in LIBRARY.rglob("*.pdf"):
        rel = p.relative_to(LIBRARY)
        if rel.parts and rel.parts[0] == "orphans":
            continue
        out[p.name].append(p)
    return out


def load_bib() -> list[dict]:
    parser = bibtexparser.bparser.BibTexParser(common_strings=True)
    with BIB.open() as f:
        db = bibtexparser.load(f, parser=parser)
    return db.entries


def get_year(entry: dict) -> int | None:
    y = entry.get("year")
    if not y:
        return None
    m = re.search(r"\d{4}", y)
    return int(m.group(0)) if m else None


def report(lines: list[str], heading: str | None = None) -> str:
    if heading:
        body = f"\n{'=' * 76}\n{heading}\n{'=' * 76}\n"
    else:
        body = ""
    return body + "\n".join(lines) + ("\n" if lines else "")


def md5_of(path: Path) -> tuple[str, str | None, str | None]:
    """Return (path-as-str, md5-hexdigest, error). Streams the file in chunks."""
    h = hashlib.md5()
    try:
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return (str(path), h.hexdigest(), None)
    except OSError as e:
        return (str(path), None, str(e))


def pdfinfo_of(path: Path) -> tuple[str, int | None, str | None]:
    """Return (path-as-str, page-count, error). Uses ``pdfinfo`` from poppler."""
    try:
        r = subprocess.run(
            ["pdfinfo", str(path)],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except FileNotFoundError:
        return (str(path), None, "pdfinfo not installed")
    except subprocess.TimeoutExpired:
        return (str(path), None, "pdfinfo timed out")
    if r.returncode != 0:
        msg = (r.stderr or r.stdout or "").strip().splitlines()
        return (str(path), None, msg[0] if msg else f"pdfinfo exit {r.returncode}")
    for line in r.stdout.splitlines():
        if line.startswith("Pages:"):
            try:
                return (str(path), int(line.split(":", 1)[1].strip()), None)
            except ValueError:
                return (str(path), None, f"unparseable Pages line: {line!r}")
    return (str(path), None, "no Pages line in pdfinfo output")


def render_histogram_png(
    decade_records: dict[int, int],
    pre_1750_records: int,
    headline: dict[str, int],
    out_path: Path,
) -> None:
    """Render the decade histogram (1750+ plus a pre-1750 bin) with a stat overlay."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    decades = sorted(decade_records.keys())
    labels: list[str] = []
    counts: list[int] = []
    if pre_1750_records:
        labels.append("<1750")
        counts.append(pre_1750_records)
    labels.extend(f"{d}s" for d in decades)
    counts.extend(decade_records[d] for d in decades)

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(labels, counts, color="#3b6ea5", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("decade")
    ax.set_ylabel("publications in library")
    ax.set_title(f"{LIBRARY_NAME} literature collection by decade")
    ax.tick_params(axis="x", rotation=45)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    overlay = (
        f"{headline['records']:,} bib records\n"
        f"{headline['pdfs']:,} PDFs\n"
        f"{headline['total_pages']:,} total pages "
        f"(median {headline['median_pages']}/PDF)"
    )
    ax.text(
        0.02, 0.97, overlay,
        transform=ax.transAxes, va="top", ha="left",
        fontsize=11, family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="#cccccc"),
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def update_readme_block(readme_path: Path, png_rel: str, headline: dict[str, int]) -> None:
    """Replace the stats sentinel block in readme.md with fresh numbers + PNG ref."""
    text = readme_path.read_text()
    pattern = re.compile(
        rf"({re.escape(README_BEGIN)})(.*?)({re.escape(README_END)})",
        re.DOTALL,
    )
    if not pattern.search(text):
        raise SystemExit(
            f"sentinel block not found in {readme_path.name}; "
            f"add a `{README_BEGIN}` … `{README_END}` block where you want the stats."
        )
    block = (
        f"\n"
        f"![Library by decade]({png_rel})\n"
        f"\n"
        f"**{headline['pdfs']:,} PDFs · {headline['records']:,} bib records · "
        f"{headline['total_pages']:,} total pages** "
        f"(mean {headline['mean_pages']:.0f} pages/PDF, median {headline['median_pages']})\n"
    )
    new_text = pattern.sub(rf"\1{block}\3", text)
    readme_path.write_text(new_text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true", help="list every offender, not just counts")
    ap.add_argument(
        "--bib",
        default=None,
        help="path to the .bib (default: the single .bib at the library root)",
    )
    ap.add_argument(
        "--scan-pdfs",
        action="store_true",
        help="also hash each PDF and run pdfinfo (slower; ~1 min for ~2000 PDFs)",
    )
    ap.add_argument(
        "--emit-readme",
        action="store_true",
        help=f"regenerate {HISTOGRAM_PNG.relative_to(REPO)} and the stats sentinel "
             f"block in {README.name}; implies --scan-pdfs",
    )
    args = ap.parse_args()
    if args.emit_readme:
        args.scan_pdfs = True

    global BIB
    BIB = _discover_bib(args.bib)

    entries = load_bib()
    lib_pdfs = collect_library_pdfs()
    paths_by_name = pdf_paths_by_basename()

    n_with_file = sum(1 for e in entries if (e.get("file") or "").strip())
    n_with_doi = sum(1 for e in entries if (e.get("doi") or "").strip())
    n_with_url = sum(1 for e in entries if (e.get("url") or "").strip())
    n = len(entries)

    out_lines: list[str] = [
        f"{BIB.name} validation",
        f"  bib: {BIB}",
        f"  library: {LIBRARY}",
        f"  records:  {n}",
        f"  pdfs:     {len(lib_pdfs)} (excluding library/orphans/)",
        "",
        f"  with file field:  {n_with_file:>5}    without: {n - n_with_file:>5}",
        f"  with doi field:   {n_with_doi:>5}    without: {n - n_with_doi:>5}",
        f"  with url field:   {n_with_url:>5}    without: {n - n_with_url:>5}",
    ]

    # ----- 1. Index bib by file basename -----
    bib_by_file: dict[str, list[dict]] = defaultdict(list)
    bib_no_file: list[dict] = []
    for e in entries:
        f = (e.get("file") or "").strip()
        if not f:
            bib_no_file.append(e)
        else:
            # Normalise to a basename, because that is what corpus itself does:
            # `BibIndex` keys on `Path(part).name.lower()` (bib/parser.py), so
            # `file = {library/D/Dunn2005.pdf}` and `file = {Dunn2005.pdf}` are
            # the same entry to a build. Comparing the raw string here made a
            # validator stricter than the product — a generated bib that
            # `corpus run` accepts reported every file as missing AND every PDF
            # as an orphan, which reads as catastrophe rather than as a
            # difference of convention. Split on `;`/`,` for the same reason
            # BibIndex does: some dialects list several files.
            for part in re.split(r"[;,]", f):
                name = Path(part.strip()).name
                if name:
                    bib_by_file[name].append(e)

    # ----- 2. Inventory gaps -----
    gaps_lines: list[str] = []

    bib_referenced = set(bib_by_file.keys())
    missing_files = sorted(name for name in bib_referenced if name not in lib_pdfs)
    orphan_pdfs = sorted(name for name in lib_pdfs if name not in bib_referenced)
    multi_referenced = {name: ents for name, ents in bib_by_file.items() if len(ents) > 1}
    duplicate_basenames = {name: paths for name, paths in paths_by_name.items() if len(paths) > 1}

    # Unlike the siphonophore library, this bib is deliberately wider than
    # the shelf: it records every paper the harvest found, and
    # only a minority of those have an obtainable open-access PDF. An
    # entry without a `file` field is the normal case here, not a gap.
    gaps_lines.append(f"bib entries with no `file` field (expected; no OA PDF): {len(bib_no_file)}")
    gaps_lines.append(f"bib `file` references that don't exist on disk:      {len(missing_files)}")
    gaps_lines.append(f"PDFs in library/ not referenced by any bib entry:    {len(orphan_pdfs)}")
    gaps_lines.append(f"PDFs referenced by more than one bib entry:          {len(multi_referenced)}")
    gaps_lines.append(f"PDF basenames that exist in multiple library/ subdirs: {len(duplicate_basenames)}")

    if args.full or len(missing_files) <= 20:
        if missing_files:
            gaps_lines.append("\n  -- bib `file` referring to missing PDFs:")
            for name in missing_files:
                ents = bib_by_file[name]
                key = ents[0].get("ID", "?")
                gaps_lines.append(f"     {name}  (cited by {key})")
    if args.full or len(orphan_pdfs) <= 30:
        if orphan_pdfs:
            gaps_lines.append("\n  -- PDFs not referenced by any bib entry:")
            for name in orphan_pdfs:
                # show subdir to disambiguate
                paths = paths_by_name[name]
                gaps_lines.append(f"     {paths[0].relative_to(LIBRARY)}")
    if multi_referenced:
        gaps_lines.append("\n  -- PDFs cited by multiple bib entries:")
        for name, ents in sorted(multi_referenced.items()):
            keys = ", ".join(e.get("ID", "?") for e in ents)
            gaps_lines.append(f"     {name}  -> {keys}")
    if duplicate_basenames:
        gaps_lines.append("\n  -- PDF basenames present in multiple subdirs:")
        for name, paths in sorted(duplicate_basenames.items()):
            gaps_lines.append(f"     {name}: {[str(p.relative_to(LIBRARY)) for p in paths]}")

    out_lines.append(report(gaps_lines, "INVENTORY GAPS"))

    # ----- 3. Metadata health -----
    health_lines: list[str] = []

    # Duplicate citation keys
    id_counts = Counter(e["ID"] for e in entries)
    dup_keys = {k: n for k, n in id_counts.items() if n > 1}

    # Duplicate DOIs (warn — sometimes intentional; e.g. plate excerpts of one work)
    doi_to_keys: dict[str, list[str]] = defaultdict(list)
    for e in entries:
        doi = (e.get("doi") or "").strip().lower()
        if doi:
            doi_to_keys[doi].append(e["ID"])
    dup_dois = {d: ks for d, ks in doi_to_keys.items() if len(ks) > 1}

    # Required-ish fields
    missing_author = [e for e in entries if not (e.get("author") or "").strip()]
    missing_year = [e for e in entries if get_year(e) is None]
    missing_title = [e for e in entries if not (e.get("title") or "").strip()]
    missing_journal = [e for e in entries if not (e.get("journal") or "").strip()]

    # No-DOI in modern era
    modern_no_doi = [
        e for e in entries
        if (yr := get_year(e)) is not None and yr >= 2000 and not (e.get("doi") or "").strip()
    ]

    # No URL anywhere (no doi, no url field)
    no_url_no_doi = [
        e for e in entries
        if not (e.get("doi") or "").strip() and not (e.get("url") or "").strip()
    ]

    health_lines.append(f"duplicate citation keys:                {len(dup_keys)}")
    health_lines.append(f"duplicate DOIs (warning, may be ok):    {len(dup_dois)}")
    health_lines.append(f"missing author:                         {len(missing_author)}")
    health_lines.append(f"missing year:                           {len(missing_year)}")
    health_lines.append(f"missing title:                          {len(missing_title)}")
    health_lines.append(f"missing journal/source:                 {len(missing_journal)}")
    health_lines.append(f"entries from 2000+ with no DOI:         {len(modern_no_doi)}")
    health_lines.append(f"entries with neither DOI nor URL:       {len(no_url_no_doi)}")

    if dup_keys:
        health_lines.append("\n  -- duplicate keys:")
        for k, n in sorted(dup_keys.items()):
            health_lines.append(f"     {k}  ({n}x)")
    if dup_dois:
        health_lines.append("\n  -- duplicate DOIs (sometimes legitimate, e.g. plate excerpts):")
        for d, ks in sorted(dup_dois.items()):
            health_lines.append(f"     {d}: {', '.join(ks)}")
    if args.full or len(missing_year) <= 10:
        if missing_year:
            health_lines.append("\n  -- entries with no year:")
            for e in missing_year:
                health_lines.append(f"     {e['ID']}: {(e.get('title') or '')[:80]}")
    if args.full or len(missing_journal) <= 10:
        if missing_journal:
            health_lines.append("\n  -- entries with no journal/source:")
            for e in missing_journal:
                health_lines.append(f"     {e['ID']}: {(e.get('title') or '')[:80]}")
    if args.full and modern_no_doi:
        health_lines.append("\n  -- 2000+ entries with no DOI:")
        for e in modern_no_doi:
            yr = get_year(e)
            health_lines.append(f"     {e['ID']} ({yr}): {(e.get('title') or '')[:80]}")

    out_lines.append(report(health_lines, "METADATA HEALTH"))

    # ----- 4. Data integrity -----
    integrity_lines: list[str] = []

    bad_dois: list[tuple[str, str]] = []
    bad_urls: list[tuple[str, str]] = []
    doi_url_mismatch: list[str] = []
    for e in entries:
        doi = (e.get("doi") or "").strip()
        url = (e.get("url") or "").strip()
        if doi and not DOI_RE.match(doi):
            bad_dois.append((e["ID"], doi))
        if url and not URL_RE.match(url):
            bad_urls.append((e["ID"], url))
        if doi and url and url.startswith("https://doi.org/") and url[len("https://doi.org/"):] != doi:
            doi_url_mismatch.append(e["ID"])

    integrity_lines.append(f"malformed DOIs:                         {len(bad_dois)}")
    integrity_lines.append(f"malformed URLs:                         {len(bad_urls)}")
    integrity_lines.append(f"doi.org URL doesn't match `doi` field:  {len(doi_url_mismatch)}")

    for label, items in (("malformed DOIs", bad_dois), ("malformed URLs", bad_urls)):
        if items and (args.full or len(items) <= 10):
            integrity_lines.append(f"\n  -- {label}:")
            for k, v in items:
                integrity_lines.append(f"     {k}: {v}")
    if doi_url_mismatch and (args.full or len(doi_url_mismatch) <= 10):
        integrity_lines.append("\n  -- DOI/URL mismatches:")
        for k in doi_url_mismatch:
            integrity_lines.append(f"     {k}")

    out_lines.append(report(integrity_lines, "DATA INTEGRITY"))

    # ----- 5. Year stats -----
    by_decade: dict[str, dict[str, int]] = defaultdict(lambda: dict(total=0, with_doi=0, with_file=0, with_url=0))
    for e in entries:
        yr = get_year(e)
        if yr is None:
            decade = "unknown"
        elif yr < 1800:
            decade = "<1800"
        else:
            decade = f"{yr // 10 * 10}s"
        bucket = by_decade[decade]
        bucket["total"] += 1
        if (e.get("doi") or "").strip():
            bucket["with_doi"] += 1
        if (e.get("file") or "").strip():
            bucket["with_file"] += 1
        if (e.get("url") or "").strip():
            bucket["with_url"] += 1

    def decade_key(d: str) -> tuple[int, str]:
        if d == "unknown":
            return (10000, d)
        if d == "<1800":
            return (-1, d)
        return (int(d[:-1]), d)

    stats_lines: list[str] = []
    stats_lines.append(f"  {'decade':<8} {'records':>8} {'with file':>10} {'with doi':>10} {'with url':>10}")
    for d in sorted(by_decade.keys(), key=decade_key):
        b = by_decade[d]
        stats_lines.append(
            f"  {d:<8} {b['total']:>8} "
            f"{b['with_file']:>10} {b['with_doi']:>10} {b['with_url']:>10}"
        )

    out_lines.append(report(stats_lines, "STATS BY DECADE"))

    # ----- 6. PDF content scan (slow; opt-in) -----
    if args.scan_pdfs:
        pdfs = [p for p in LIBRARY.rglob("*.pdf")
                if not (p.relative_to(LIBRARY).parts
                        and p.relative_to(LIBRARY).parts[0] == "orphans")]

        # Hash + page-count in parallel.
        hash_results: dict[str, str] = {}
        hash_errors: list[tuple[str, str]] = []
        page_results: dict[str, int] = {}
        page_errors: list[tuple[str, str]] = []

        with ProcessPoolExecutor() as ex:
            future_kind: dict = {}
            for p in pdfs:
                future_kind[ex.submit(md5_of, p)] = "hash"
                future_kind[ex.submit(pdfinfo_of, p)] = "page"
            for fut in as_completed(future_kind):
                kind = future_kind[fut]
                path, value, err = fut.result()
                if kind == "hash":
                    if err is not None:
                        hash_errors.append((path, err))
                    else:
                        hash_results[path] = value  # type: ignore[assignment]
                else:  # page
                    if err is not None:
                        page_errors.append((path, err))
                    else:
                        page_results[path] = value  # type: ignore[assignment]

        # Group by hash to find content duplicates.
        by_hash: dict[str, list[str]] = defaultdict(list)
        for path, h in hash_results.items():
            by_hash[h].append(path)
        dup_hash_groups = {h: ps for h, ps in by_hash.items() if len(ps) > 1}

        scan_lines: list[str] = []
        scan_lines.append(f"PDFs scanned:                           {len(pdfs)}")
        scan_lines.append(f"successfully hashed:                    {len(hash_results)}")
        scan_lines.append(f"duplicate-content groups (same hash):   {len(dup_hash_groups)}")
        scan_lines.append(f"corrupt / unreadable (pdfinfo failed):  {len(page_errors)}")
        scan_lines.append(f"hash read errors:                       {len(hash_errors)}")

        if page_results:
            pages = list(page_results.values())
            scan_lines.append("")
            scan_lines.append(f"total pages:                            {sum(pages)}")
            scan_lines.append(f"mean pages per pdf:                     {statistics.mean(pages):.1f}")
            scan_lines.append(f"median pages per pdf:                   {statistics.median(pages):.0f}")

        if dup_hash_groups:
            scan_lines.append("\n  -- duplicate-content PDF groups:")
            for h, paths in sorted(dup_hash_groups.items(), key=lambda kv: kv[1][0]):
                scan_lines.append(f"     {h[:12]}…  ({len(paths)} files)")
                for path in sorted(paths):
                    rel = Path(path).relative_to(LIBRARY)
                    scan_lines.append(f"        {rel}")

        if page_errors and (args.full or len(page_errors) <= 30):
            scan_lines.append("\n  -- pdfinfo failures (corrupt or unreadable):")
            for path, err in sorted(page_errors):
                rel = Path(path).relative_to(LIBRARY)
                scan_lines.append(f"     {rel}: {err}")

        if hash_errors and (args.full or len(hash_errors) <= 30):
            scan_lines.append("\n  -- hash read errors:")
            for path, err in sorted(hash_errors):
                rel = Path(path).relative_to(LIBRARY)
                scan_lines.append(f"     {rel}: {err}")

        out_lines.append(report(scan_lines, "PDF CONTENT SCAN"))

    # ----- 7. Emit README assets (PNG histogram + sentinel block) -----
    if args.emit_readme:
        # Bin records by decade for the histogram. Pre-1750 entries are summarized
        # in the overlay rather than shown on the axis.
        decade_records: dict[int, int] = defaultdict(int)
        pre_1750 = 0
        for e in entries:
            yr = get_year(e)
            if yr is None:
                continue
            if yr < 1750:
                pre_1750 += 1
                continue
            decade_records[yr // 10 * 10] += 1

        pages = list(page_results.values()) if page_results else []
        headline = {
            "records": n,
            "pdfs": len(lib_pdfs),
            "total_pages": sum(pages),
            "mean_pages": statistics.mean(pages) if pages else 0,
            "median_pages": int(statistics.median(pages)) if pages else 0,
        }

        ASSETS.mkdir(exist_ok=True)
        render_histogram_png(decade_records, pre_1750, headline, HISTOGRAM_PNG)
        png_rel = HISTOGRAM_PNG.relative_to(REPO).as_posix()
        update_readme_block(README, png_rel, headline)
        print(f"\nwrote {HISTOGRAM_PNG.relative_to(REPO)}", file=sys.stderr)
        print(f"updated stats block in {README.name}", file=sys.stderr)

    output = "\n".join(out_lines) + "\n"
    OUT.write_text(output)
    sys.stdout.write(output)


if __name__ == "__main__":
    main()

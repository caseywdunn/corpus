#!/usr/bin/env python3
"""Grow the bib from what the library already holds — mine the reference lists
of papers you have for papers you do not.

TEMPLATE — copied into a library repo by the `assemble-library` skill.

This is the step that makes a library deeper than an index query. A 19th-century
paper usually enters a bibliography only because someone later cited it, so the
reference lists of papers you hold are the only route to the literature no index
will surface.

## Read Grobid's output. Do not re-parse the PDFs.

**This runs after a first `corpus run`, not during assembly.** That ordering is
the whole design, and it is worth understanding before changing it.

`corpus run` already sends every PDF through Grobid, which is a purpose-built
reference parser, and leaves the result at
`<output_dir>/documents/<hash>/grobid.tei.xml`. Re-extracting references from
the PDF text with a regex duplicates that work and does it far worse. Measured
on the four demo papers, Grobid against a regex over `pdftotext`:

    grobid= 15   regex=  2
    grobid= 15   regex=  7
    grobid= 26   regex=  3
    grobid=  1   regex=  0

Roughly 20% recall, and what the regex does return are fragments — continuation
lines like "Eur J Neurosci 18: 1848-1860" rather than whole references. Reference
lists wrap across lines and columns in layouts that vary by publisher and by
century; a regex loses to a trained model here and always will.

So the pipeline for a library is: assemble, build once, mine, add, build again.
That is a lap of the corpuscle-improvement loop, and mining is one of the things
a lap is *for*.

## The gate goes on the resolved record, not the reference string

A clade paper's reference list is mostly general biology, statistics and
methods. Testing the raw reference string for clade vocabulary lets through
anything whose title happens to mention the group and rejects everything whose
relevance is only visible once resolved. Resolve first, then judge.

## Where the vocabulary comes from

The clade's own taxonomy. `taxonomy.dwca.zip` lists every genus and species in
the group, including the synonyms and segregate genera the older literature
used — a better vocabulary than anything hand-written, and it stays correct when
the snapshot is rebuilt. `relevance.yaml` adds what a taxonomy cannot know:
vernacular names, former family placements, subject terms.

## Beyond this script

Once a corpuscle is *served*, `get_missing_references` answers the same question
with better evidence: the citation graph reconciled corpus-wide, deduplicated,
and ranked by how often each absent work is cited. This script is the local,
pre-serve version. Prefer the served one when you have it.

Usage:

    python mine_references.py --output-dir ../output
    python mine_references.py --output-dir ../output --min-citations 2
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

try:
    import yaml
except ImportError:
    sys.exit("pyyaml is required: it is in the library environment.yaml")

REPO = Path(__file__).resolve().parents[1]
BUILD = REPO / "build"
OUT = BUILD / "mined_references.jsonl"

TEI_NS = {"tei": "http://www.tei-c.org/ns/1.0"}
_DOI = re.compile(r"\b10\.\d{4,9}/\S+", re.I)


def clade_vocabulary(taxonomy_zip: Path, extra: Path) -> set[str]:
    """Lowercased terms that make a resolved record plausibly on-topic.

    Genus-level and above only: species epithets alone ("elegans", "vulgaris")
    are ordinary words and would let almost anything through.
    """
    vocab: set[str] = set()
    if taxonomy_zip.exists():
        with zipfile.ZipFile(taxonomy_zip) as zf:
            names = [n for n in zf.namelist()
                     if n.lower().endswith(("taxon.tsv", "taxa.tsv"))]
            if names:
                with zf.open(names[0]) as fh:
                    reader = csv.DictReader(
                        io.TextIOWrapper(fh, "utf-8", errors="replace"), delimiter="\t"
                    )
                    for row in reader:
                        name = (row.get("scientificName") or "").strip()
                        rank = (row.get("taxonRank") or "").strip().lower()
                        if not name:
                            continue
                        first = name.split()[0]
                        if len(first) > 3 and first[0].isupper():
                            vocab.add(first.lower())
                        if rank in {"genus", "family", "order", "class", "phylum"}:
                            vocab.add(name.lower())
    if extra.exists():
        spec = yaml.safe_load(extra.read_text(encoding="utf-8")) or {}
        for term in spec.get("terms") or []:
            if str(term).strip():
                vocab.add(str(term).strip().lower())
    return vocab


def is_relevant(record: dict, vocab: set[str]) -> bool:
    """Does this *resolved* record belong in the library?

    Generous about where a term appears — a paper can be about the clade without
    naming it in the title — and strict about requiring one at all.
    """
    if not vocab:
        return True  # nothing configured: let a human filter instead
    haystack = " ".join(
        str(x).lower() for x in (
            record.get("title") or "",
            record.get("journal") or "",
            " ".join(record.get("authors") or []),
        )
    )
    return any(term in haystack for term in vocab)


def references_from_tei(tei_path: Path) -> list[dict]:
    """Every parsed reference in one document's Grobid TEI.

    Grobid has already done the hard part — segmenting the reference list and
    parsing each entry into fields. This reads that, and does not second-guess
    it.
    """
    try:
        root = ET.parse(tei_path).getroot()
    except (ET.ParseError, OSError):
        return []

    out: list[dict] = []
    back = root.find(".//tei:back", TEI_NS)
    if back is None:
        return out

    for bibl in back.iterfind(".//tei:biblStruct", TEI_NS):
        title_el = bibl.find(".//tei:title[@level='a']", TEI_NS)
        if title_el is None:
            title_el = bibl.find(".//tei:title", TEI_NS)
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""

        journal_el = bibl.find(".//tei:title[@level='j']", TEI_NS)
        journal = "".join(journal_el.itertext()).strip() if journal_el is not None else ""

        authors = []
        for pers in bibl.iterfind(".//tei:author/tei:persName", TEI_NS):
            surname = pers.find("tei:surname", TEI_NS)
            if surname is not None and (surname.text or "").strip():
                authors.append(surname.text.strip())

        year = None
        date = bibl.find(".//tei:date[@when]", TEI_NS)
        if date is not None:
            m = re.search(r"\b(1[5-9]\d{2}|20[0-2]\d)\b", date.get("when", ""))
            if m:
                year = int(m.group(1))

        doi = None
        idno = bibl.find(".//tei:idno[@type='DOI']", TEI_NS)
        if idno is not None and (idno.text or "").strip():
            doi = idno.text.strip()
        else:
            m = _DOI.search("".join(bibl.itertext()))
            if m:
                doi = m.group(0).rstrip(".,;)")

        if not title and not doi:
            continue  # nothing to search for; a mis-parsed entry, not a lead
        out.append({
            "title": title, "journal": journal, "authors": authors,
            "year": year, "doi": doi,
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", default=None,
                    help="built corpuscle output dir (default: ../output)")
    ap.add_argument("--taxonomy", default=None, help="default: taxonomy.dwca.zip")
    ap.add_argument("--min-citations", type=int, default=1,
                    help="only report references cited by at least N held papers")
    ap.add_argument("--all", action="store_true",
                    help="skip the relevance gate (report everything cited)")
    args = ap.parse_args()

    outdir = Path(args.output_dir) if args.output_dir else REPO / "output"
    docs = outdir / "documents"
    if not docs.is_dir():
        sys.exit(
            f"no built corpuscle at {outdir}.\n"
            "This step reads Grobid's parsed reference lists, which `corpus run`\n"
            "produces — so build once first, then mine. Re-parsing the PDFs here\n"
            "instead recovers about a fifth of the references, in fragments."
        )

    tax = Path(args.taxonomy) if args.taxonomy else REPO / "taxonomy.dwca.zip"
    relevance = REPO / "relevance.yaml"
    vocab = clade_vocabulary(tax, relevance)
    print(f"clade vocabulary: {len(vocab):,} terms", file=sys.stderr)
    if not vocab and not args.all:
        print("  ! empty — every reference will pass the relevance gate.",
              file=sys.stderr)
    elif not relevance.exists():
        # Measured on the demo set: taxonomy names alone rejected core
        # literature because the papers print "Siphonophoren" where the
        # taxonomy says "Siphonophorae". Warn before the damage, not after.
        print(
            "  ! no relevance.yaml — the vocabulary is formal taxon names only.\n"
            "    The literature prints vernacular and inflected forms, so this\n"
            "    gate will reject on-topic papers. See relevance.yaml.example.",
            file=sys.stderr,
        )

    teis = sorted(docs.glob("*/grobid.tei.xml"))
    if not teis:
        sys.exit(
            f"no grobid.tei.xml under {docs}. Was Grobid running during the "
            "build? Without it there are no parsed reference lists to mine."
        )

    by_key: dict[str, dict] = {}
    citations: Counter = Counter()
    n_refs = 0
    for tei in teis:
        for ref in references_from_tei(tei):
            n_refs += 1
            key = (ref["doi"] or ref["title"][:120]).lower()
            citations[key] += 1
            by_key.setdefault(key, ref)

    kept = [
        {**ref, "cited_by": citations[key]}
        for key, ref in by_key.items()
        if citations[key] >= args.min_citations
        and (args.all or is_relevant(ref, vocab))
    ]
    kept.sort(key=lambda r: (-r["cited_by"], r.get("title") or ""))

    BUILD.mkdir(exist_ok=True)
    with OUT.open("w", encoding="utf-8") as sink:
        for ref in kept:
            sink.write(json.dumps(ref, ensure_ascii=False) + "\n")

    print(
        f"\n{n_refs:,} references across {len(teis):,} documents, "
        f"{len(by_key):,} distinct, {len(kept):,} kept "
        f"-> {OUT.relative_to(REPO)}",
        file=sys.stderr,
    )
    for ref in kept[:5]:
        print(f"    {ref['cited_by']}x  {(ref.get('title') or '')[:72]}", file=sys.stderr)

    # An over-tight vocabulary is the quiet failure here: it does not error, it
    # just silently discards the literature you were looking for. Surface it.
    if by_key and not args.all:
        rejected = len(by_key) - len(kept)
        if rejected / len(by_key) > 0.4:
            print(
                f"\n  ! the relevance gate rejected {rejected}/{len(by_key)} "
                f"({rejected / len(by_key):.0%}).\n"
                "    Re-run with --all and read a sample of what was dropped. A "
                "vocabulary\n    missing the inflections the literature uses "
                "looks exactly like this.",
                file=sys.stderr,
            )
    print(
        "\nThese are leads, not entries. Resolve each against Crossref before "
        "adding it to the bib — a Grobid-parsed reference is as good as the "
        "citing author's, and no better.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

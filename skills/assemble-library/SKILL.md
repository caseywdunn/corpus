---
name: assemble-library
description: >
  Assemble a corpus-ready library of scientific literature for a clade or
  taxonomic group — harvest the bibliography from public indexes, retrieve
  the open-access PDFs, build a Darwin Core taxonomy snapshot, and write the
  config.yaml that `corpus run` needs. Use this whenever someone wants to
  gather, harvest, or build a literature collection for a taxon, start a new
  corpus or corpuscle for a group of organisms, or turn "I want a corpus of
  everything written about X" into files on disk — even if they don't say
  "corpus" or "library". Requires the `corpus` CLI installed; does NOT
  require a running MCP server.
---

# Assemble a library for a clade

A **library** is the upstream collection: PDFs, a `.bib`, a taxonomy snapshot
and the clade-specific knowledge files. A **corpuscle** is what `corpus run`
builds from it. This skill produces the first so the second becomes one command.

The output is a directory `corpus run` accepts with no further setup:

```text
<clade>/
├── library/              # PDFs, sharded by surname-letter shelf
│   └── orphans/          # PDFs kept but with no known bibliographic record
├── <clade>.bib           # every paper discovered; file = {...} only where a PDF exists
├── lexicon.yaml          # DRAFT — needs expert review
├── instructions.md       # DRAFT — needs expert review
├── taxonomy.dwca.zip
├── config.yaml
├── readme.md
├── CONTRIBUTING.md       # the pipeline order for this library, generated
├── environment.yaml
├── .gitattributes        # Git LFS for *.pdf
└── scripts/              # the harvest/curation pipeline, copied and editable
```

## The one idea that shapes everything else

**The bib is a superset of the PDFs.** Every relevant paper you discover gets an
entry whether or not you obtained a file. Entries backed by a PDF carry
`file = {...}`; everything else is a **want-list** that stays queryable — corpus
will answer "what are the most-cited papers I don't have?" from it, and it is
what lets someone add PDFs later.

This matters because retrieval mostly fails, and that is normal. Expect roughly
a quarter to a third of fetch attempts to succeed. A failure is a want-list
entry, not a defeat. Say this to the user early so a 30% hit rate reads as the
expected outcome rather than a broken run.

## Keys and contact address — settle these at the start

Two environment variables decide how much of the literature is reachable, and
both are cheaper to set now than to discover missing halfway through a harvest.

| Variable | Effect if unset |
|---|---|
| `CORPUS_CONTACT_EMAIL` | Slower harvest, more 429s. Not a secret; always set it |
| `BHL_API_KEY` | **The pre-1930 literature does not happen** |

That second one deserves saying plainly rather than listing: BHL is the only
practical route to article-level records for older material, so for a clade with
a deep historical literature an absent key is the difference between a corpus
starting in 1990 and one starting in 1758. Ask for it before harvesting, point
the user at <https://www.biodiversitylibrary.org/getapikey.aspx> (free, arrives
by email), and if they would rather not wait, say what will be missing and
proceed.

Never accept a key pasted into the conversation — it lands in the transcript.
`references/retrieval-ethics.md` has the safe ways to set one, and the rules for
handling keys generally.

Report at the end which keys were absent. A skipped source that nobody mentions
is indistinguishable from a literature that does not exist.

## Before anything else: scope the clade with the user

**Do not skip this and do not guess.** The inclusion rule is the single
highest-leverage decision in the build, it is genuinely contentious, and it is
much cheaper to agree now than to re-harvest later.

Settle five things, and write the agreed rule into the generated
`build_bib.py` docstring so it stays discoverable:

1. **The inclusion rule.** What makes a paper "about" this group? Viburnum's is
   *"genus or vernacular name in title or abstract; family- and order-level
   systematics; Caprifoliaceae treatments up to 2003 but not after, because APG
   II moved the genus out of that family."* Note how much clade-specific
   judgment is in that one sentence.
2. **Taxonomic authority.** Which checklist is the truth for this group? See
   `references/taxonomy-sources.md`.
3. **Historical and segregate names.** What did the older literature publish
   under? A 19th-century paper on *Viburnum tinus* may only ever write
   *Tinus laurifolius*, and no modern-name query will find it.
4. **Vernacular names**, including non-English ones — often the only way into
   regional literature.
5. **Homonym traps.** Ask for these explicitly and propose candidates for the
   user to confirm, because they are subtle and they bite. "Kalina" is Slavic
   for *Viburnum opulus* and also a thermodynamic power cycle with thousands of
   engineering papers. "Viburnum Trend" is a Missouri lead-mining district. A
   vernacular sweep pulls in both.

When a query turns out to be contaminated, quarantine it rather than deleting
the cache: accept its records only where the title is independently on-topic.

## Retrieval ethics — non-negotiable

Read `references/retrieval-ethics.md` before writing any fetch code, and carry
its substance into the generated library's own docs.

The short form: use documented public APIs and the OA endpoints publishers
advertise. Honour `robots.txt`, rate limits and `Retry-After`. Identify the
client with a real contact address, read from `CORPUS_CONTACT_EMAIL` and never
hardcoded.

Ask the user to set that before the harvest begins, and say why — Crossref,
Unpaywall and OpenAlex all give identified clients materially higher rate
limits, so an unset variable means a slower harvest and more 429s:

```bash
export CORPUS_CONTACT_EMAIL="you@example.edu"
```

Scripts should fail with that instruction rather than fall back to a default;
there is no safe default, since a placeholder lies to the service and someone
else's address sends them your traffic.
Do not defeat bot protection, rotate or spoof a User-Agent to evade a block, use
shadow-library mirrors, or use anyone's institutional credentials.

**A 403 is a final answer, not a puzzle.** Record it and move on. This is not
hypothetical: BHL's `partpdf` and `itempdf` endpoints refuse non-interactive
clients whatever you send, so those entries land in the bib with an id, a
`rights` statement and no `file`. Publisher 403s are the dominant failure mode
and are simply not engineered around.

## The pipeline

Steps 3–5 are a loop. Mined references become bib entries, entries become
download targets, downloaded PDFs yield more references. **Two passes is where
it stops paying.**

```bash
python scripts/openalex_probe.py --check                # 0. preflight
python scripts/harvest_openalex.py --dry-run            #    cost it first

corpus taxonomy ingest --source dwca \
    --input <checklist.zip> --root-id <taxonID>         # 1. taxonomy
corpus taxonomy export -o taxonomy.dwca.zip

python scripts/harvest_openalex.py                      # 2. harvest
python scripts/harvest_bhl.py                           #    needs BHL_API_KEY
python scripts/harvest_bhl.py --fulltext                #    for historical names

python scripts/build_bib.py --report                    # 3. assemble

python scripts/resolve_oa.py                            # 4. locate + fetch
python scripts/fetch_pdfs.py

python scripts/mine_references.py --source both         # 5. grow
python scripts/build_bib.py
python scripts/fetch_pdfs.py

python scripts/validate_bib.py --emit-readme            # 6. check
```

Intermediate state belongs in `build/`, gitignored and re-derivable — but
`build/records.jsonl` costs a day's OpenAlex quota to rebuild, so do not delete
it casually.

### 1. Taxonomy

Build the Darwin Core snapshot first; everything downstream resolves names
against it.

**This needs no script — `corpus` already does it.** Two paths, and both end in
a committed `taxonomy.dwca.zip` so later runs ingest locally and offline rather
than re-walking a rate-limited API:

```bash
# Marine clade — walk WoRMS from an AphiaID
corpus taxonomy ingest --source worms --root-id <AphiaID>

# Everything else — download a checklist archive, prune it to your clade
corpus taxonomy ingest --source dwca --input <checklist.zip> --root-id <taxonID>

corpus taxonomy export -o taxonomy.dwca.zip     # either path
```

`--root-id` prunes a full checklist to the subtree you care about, so you can
download a large archive and keep only your clade — verified at 801 records in,
66 out. That is why no bespoke fetcher is needed: the existing library's
`build_taxonomy.py` walks GBIF's checklist API to avoid downloading a whole
archive, and download-then-prune reaches the same place with a tool that already
exists and is tested.

**Prefer one curated checklist over the GBIF backbone**, which merges sources
and leaves unplaced combinations dangling with no synonymy. Per-source gotchas,
and what does *not* survive ingest, are in `references/taxonomy-sources.md`.

### 2. Harvest — broad, and before you filter

Run several independent index sweeps and **tag every record with the query that
found it**. The tag is evidence: a record found by a title search is a different
kind of hit from one found by an outgroup query, and the relevance filter should
be able to reason over *how* something was found.

Sources, and what each is for, are in `references/harvesting.md`. In short:
OpenAlex for modern literature, BHL Part records for pre-1930 material searched
under historical and segregate names, Crossref / PubMed / Europe PMC as the
clade warrants, plus any lab or author archive the user names.

The OpenAlex queries live in `scripts/queries.yaml`, **not** in the harvest
script — they are the output of the scoping conversation, so someone revising the
harvest for a new clade edits a list of queries rather than a program. Copy
`queries.yaml.example` and rewrite it. `harvest_openalex.py` is resumable per
query and per page, deduplicates on the work id, and keeps the tag of whichever
query found a record first.

**Budget the OpenAlex harvest before running it.** It meters by credit, and the
free allowance is about 100 search requests per day — small enough that a
careless harvest spends days rather than minutes. `openalex_probe.py --estimate
'<filter>'` costs one request and tells you what the whole harvest would cost,
because `meta.count` is returned on the first page. Always page at
`per_page=200`: cost is per request, so a smaller page size multiplies the bill
for identical results. The harvest must be incremental and resumable regardless,
since a 429 here has a `Retry-After` measured in hours.

### 3. Assemble the bib

Apply the inclusion rule from the scoping conversation. Keys and filenames
follow `Smith1998` / `SmithJones1998` / `Smithetal1998`, with `Smith_Jones1998.pdf`
for the file.

**Read `references/bib-conventions.md` before implementing key generation.**
There is a specific trap there that has already caused silent data corruption in
an existing library: positional `a`/`b`/`c` suffixes shift when a duplicate
merges, which can rebind an existing PDF onto a different paper. Derive stems
from something stable, or keep a committed lockfile mapping filename → work.

### 4. Retrieve

Resolve OA locations, then try each in turn. Keep the first response that
actually starts with `%PDF`, clears a size floor, and parses under `pdfinfo` —
publisher interstitials return HTTP 200 with HTML bodies and will otherwise land
in `library/` looking like papers.

Cache failures so a re-run does not re-hammer dead ends, and make retry explicit
(`--skip-known-failures`) for when a resolver has genuinely improved.

### 5. Mine references

This is what makes the library deeper than an index query, and it is why
retrieval comes first. Two channels: Crossref reference deposits for entries
with a DOI, and `pdftotext` over `library/` for everything else — which is where
the pre-DOI literature lives, since a 19th-century paper usually enters a
bibliography only by way of someone citing it.

Apply the relevance test **to the resolved record, not the reference string**. A
clade paper's reference list is mostly general biology, statistics and methods.

### 6. Validate

`validate_bib.py` reports six sections; a clean run has `INVENTORY GAPS` empty.
Verify entries against external authorities where you can — Crossref by DOI or
bibliographic match, BHL for historical material, OpenAlex ids — and **record
how each was verified**, because provenance differs enormously. Entries
reconstructed from a reference string in a paper you happen to hold have titles
as good as the citing author's and no better; the bib should say so rather than
presenting them as equally solid.

## Generate the clade-knowledge files as review-flagged drafts

`lexicon.yaml` and `instructions.md` capture what the corpus needs to know about
this group: morphology and recurring subjects for the lexicon; nomenclatural
history, what is *not* in the clade despite older literature lumping it in, and
name collisions a search will surface, for the instructions.

Both carry a visible header saying they are agent drafts pending expert review,
and the generated `readme.md` lists reviewing them as an open task. You will get
these thin — that is expected and fine. A draft the user corrects beats an empty
file, and corpus tolerates both being imperfect.

**Write a lexicon the matcher can actually use.** Matching is exact, whole-word
and case-insensitive — not stemmed. English plurals and variants must be listed
under `synonyms`, because generated English inflection was measured and rejected
(suffixing `cnida` matched `Cnidaria` 4,444 times). Non-English forms go under
`translations`, where an inflection table does reach them. A draft of bare
canonical terms with no synonyms will match far less than it appears to.

Whether the draft is any *good* is measured after a build, not here.

## Finish

Wire up `config.yaml`, confirm with `corpus check`, and report:

- the want-list size and where to find it — this is the headline number, and it
  is a feature;
- which API keys were set and what each unlocked, so the user knows what a
  second run with more keys would add;
- that `lexicon.yaml` and `instructions.md` need expert review.

Then tell them the next step: `corpus run` inside this directory builds the
corpuscle, and the library directory *is* the corpuscle directory.

## What ships in `scripts/`

`validate_bib.py` and `environment.yaml` are **templates copied into the library
repo** — clade-agnostic, tested, and yours to edit once copied. `validate_bib.py`
discovers the `.bib` rather than naming it, so nothing needs editing for a new
group.

The library gets **its own conda environment**, separate from corpus's. It needs
`bibtexparser`, which corpus does not use (it has its own parser), and it does
none of the OCR or embedding work that makes the corpus environment heavy. So a
curator can run the harvest without a CUDA-capable torch. Activate the library
env for these scripts; activate corpus's for `corpus check` and `corpus run`.

The harvest scripts are **written per clade, not copied** — `build_bib.py` in
particular, because the inclusion rule from the scoping conversation lives in its
docstring and its tallies. Use the pipeline above as the shape and
`references/harvesting.md` for what each source is for.

## Reference files

- `references/retrieval-ethics.md` — the full rules, and why each exists. Read before writing fetch code.
- `references/taxonomy-sources.md` — choosing and building the DwC-A; per-source gotchas.
- `references/harvesting.md` — index sweeps, query tagging, rate limits, homonym quarantine.
- `references/bib-conventions.md` — keys, filenames, shelves, and the suffix-rebinding trap.

# Choosing and building the taxonomy snapshot

The Darwin Core Archive is built **first**, because every later step resolves
names against it — the relevance filter, the reference miner, and eventually
`corpus run` itself.

Export to a committed `taxonomy.dwca.zip` rather than fetching at build time.
Later runs then ingest from a local file in seconds instead of re-walking a
rate-limited REST API, and the snapshot is versioned with the library.

## Which source

| Clade | Source | How |
|---|---|---|
| Marine groups | **WoRMS** | `corpus taxonomy export` with the AphiaID subtree root |
| Vascular plants | **WCVP** (Kew) | via GBIF's checklist API; this is what POWO publishes |
| Everything else | a curated DwC-A | GBIF, ITIS or Catalogue of Life |

**Prefer one curated checklist over the GBIF backbone.** The backbone merges
many checklists, and the merge is where synonymy goes to die: for *Viburnum* it
leaves unplaced *Opulus* and *Tinus* combinations sitting under the genus with
no resolution at all. A single curated treatment is opinionated, and that is the
point — you want someone's considered answer, not the union of everyone's.

## Two details that are easy to get wrong

**Take names from `scientificName` minus the authorship, not from
`canonicalName`.** Only the former keeps the `var.` / `subsp.` connector that
the literature actually prints. Strip authorship yourself rather than accepting
a pre-canonicalised string that has quietly dropped rank connectors.

**Fetch the ancestors above your root separately and append them.** Otherwise
the root row's `parentNameUsageID` dangles outside the archive, and a consumer
walking the parent chain hits a missing id rather than a clean stop.

## What survives ingest

`corpus taxonomy ingest` keeps 15 DwC terms. Several fields a taxonomist would
expect are **not** among them — `nomenclaturalStatus`, `originalNameUsageID`,
`namePublishedIn`, `taxonRemarks`. The practical consequence: the snapshot
records *that* a name is unaccepted but not *why*.

Two things follow. Do not plan a workflow that reasons about nomenclatural
status from the snapshot alone. And when the source offers those fields, note in
the library's readme that they were available and dropped, so a later widening of
the ingest is a re-run rather than a re-investigation.

Vernacular names are a related gap: the exporter writes a vernacular extension
but no ingest path produces vernacular rows, so that extension is currently
always empty on a corpus-built snapshot. Do not rely on it.

## Authorship conventions differ by kingdom

Zoological authorship carries a year (`Linnaeus, 1758`); botanical (ICN) does
not (`Rehder`, `(Kache) Hesse`, `(Vent.) P.Silva`). The corpus authority build
detects which convention a snapshot uses and records the verdict, so downstream
tools can say "unsupported for this taxonomy" rather than silently returning
nothing.

Worth checking on a fresh snapshot: some sources map a year column into the
authorship field, which can make a botanical checklist look zoological. If the
convention verdict surprises you, inspect a few authorship strings directly
before trusting it.

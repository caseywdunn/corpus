# Harvesting — sources, tagging, and what goes wrong

Harvest **broadly and before you filter.** Run several independent index sweeps
and keep everything; relevance is decided later, in `build_bib.py`, where the
rule is written down and the tallies are inspectable.

## Tag every record with the query that found it

This is the part most likely to be skipped and most worth keeping. The tag is
*evidence about how a record was found*, and the relevance filter should reason
over it: a hit from `title.search:<genus>` is a different kind of claim from a
hit from an outgroup query or a vernacular sweep.

It is also what makes a contaminated query recoverable. When you discover a
sweep pulled in junk, you can accept records carrying that tag **only** where the
title is independently on-topic — a quarantine — instead of rebuilding the whole
cache. If the cache is ever rebuilt from scratch with the query narrowed, the
quarantine becomes a no-op and can go.

## Sources

**OpenAlex** — the modern literature. Cache every hit, deduplicated on the work
id.

> **Hard daily credit limit.** A full harvest of every query costs most of the
> free allowance; when it runs out you get 429s with a `Retry-After` measured in
> *hours*, not seconds. So the harvest must be incremental and resumable —
> re-running the next day fills in what the previous run missed. This is also why
> the record cache is precious: deleting it costs a day, not a minute.

**BHL** — where the pre-1930 material is. Collect **Part** records (article-level
segments with their own title, authors, page range and often a minted DOI) rather
than whole volumes: a 600-page flora that mentions your genus once is not a paper
about your genus. Search under the segregate and historical genera as well as the
current name, since that is what the older literature published under. Needs
`BHL_API_KEY`.

**Crossref / PubMed / Europe PMC** — as the clade warrants.

**Lab or author archives the user names.** Often the highest-yield source for a
well-studied group, because a lab that posts its own PDFs supplies material that
is paywalled everywhere else. Prefer structured metadata where the page offers it
(COinS, embedded citations) over scraping rendered text.

## Homonym traps

Ask about these during scoping, propose candidates, and have the user confirm
before harvesting. They are subtle, they are specific to the clade, and a
vernacular sweep is where they enter.

Two real ones from the viburnum library:

- **"kalina"** is Slavic for *Viburnum opulus*. It is also the name of a
  thermodynamic power cycle with a few thousand engineering papers.
- **"Viburnum Trend"** is a lead-mining district in Missouri.

The failure is not that these appear — it is that they appear in the *hundreds*
and look like a productive query. Check the top titles of any vernacular sweep by
eye before trusting its yield.

## Excluding what is not literature

Decide these with the user, because they are clade-specific and they are large:

- **`dataset` records** — for *Viburnum* these were ~7,800 GBIF occurrence
  records and GenBank accessions. Not literature.
- **Nursery and seed catalogues** — they reach a harvest through BHL because they
  list stock. Trade ephemera.
- **Marginal abstract-only mentions** — a paper naming the taxon only in a species
  list is thin, but the alternative (title-only matching) loses most of the
  community-ecology literature. This is the main source of marginal entries and
  there is no clean answer; make the call explicitly and record it.

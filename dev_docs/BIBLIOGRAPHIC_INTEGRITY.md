# Bibliographic identity and source authority

Bibliographic materialization belongs to the build plane. The served authority
DB is read-only; fixing a citation means correcting its source or build rules,
then rebuilding and bundling the results.

## Canonical fields and provenance

A library BibTeX entry supplied during metadata extraction has the same
bibliographic authority as an explicit `corpus bib import`. The authority DB
retains its entry key and supplied fields in `work_bib_sources`; its document
membership retains `extraction_method: bib` and `bib_key`. Reconciliation moves
these sources with the work and materializes the complete ordered author list.
It never unions authors by surname. Raw extracted reference observations and
reconciliation decisions remain independently auditable.

An explicit import supersedes the build-time copy with the same BibTeX key.
When independent curated entries disagree, explicit imports precede build
sources and source IDs break ties deterministically. The citation payload's
`bibliographic_conflicts` names the field, selected source, policy, and all
conflicting supplied values. This is a review signal; resolving it belongs in
the library. Extracted fields cannot win over a supplied curated value.

## Publication locators

`volume`, `number` (issue), `pages`, `eid`, `articleno`, `chapter`, `booktitle`,
`publisher`, `edition`, and `series` pass through metadata, authority storage,
BibTeX export/import, and the citation's structured `fields`. Absent values
remain null; the source's spelling and ranges are retained. Physical PDF
selection is `keeppages` and never substitutes for publication `pages`.

The current `author-year` style renders the journal or book title, volume and
issue, chapter, page range, article identifier, and publisher when supplied.
Only the rendered page range converts BibTeX `--` to an en dash. `eid` takes
precedence over `articleno` in the display when both exist; both remain in the
structured fields. Edition and series are available as structured fields but
are not currently displayed by this style.

## Updating existing builds

Normal `corpus run` invalidates metadata produced before the
`bibliographic-metadata-v2` receipt and regenerates BibTeX metadata, including
locators. The next authority pass migrates the SQLite schema and re-ingests
metadata whose consumed contents changed. Reference mappings are rematerialized
when their producer or corpus identity changes. Explicit imports continue to
be retained inputs; make the same curation changes in the library `.bib` so a
clean build reproduces them. Rebuild the served bundle after updating the build.
A server-only upgrade cannot restore fields missing from an old bundle.

## Citation units and shared identifiers

A DOI may describe a whole book, while the library holds separate volumes or
an atlas. Phase 1 compares the full available title, journal, year, edition and
publication locators before grouping documents. Without a DOI it also compares
the ordered surname list. Compatible duplicate scans still share one work.
Conflicting parts receive deterministic `#part:` work-ID suffixes; their DOI
remains intact and `shared_identifier` records the common DOI or short key.
A reference that supplies enough metadata selects the corresponding part. A
shared DOI alone stays a separate unresolved part-level match, rather than
acquiring one volume's identity arbitrarily. The document hash continues to
select exactly the requested PDF, including its local rights directives.

Reconciliation checks curated title-to-title evidence before considering
first-page text or citation popularity. Conflicting titles, explicit parts,
years, DOIs and locators cannot replace a BibTeX-backed identity. Rejected
candidates and identity migrations are recorded with their evidence in
`work_identity_decisions`. This guard still allows the existing reconstruction
of uncurated misparsed headers and strongly supported duplicate titles.

On the first authority pass using `document-identity-v2`, memberships are
re-derived from current metadata and the reference graph is rematerialized.
Previously conflated documents may acquire new work IDs; callers should resolve
saved document hashes again. `work_reconciliation_decisions` and immutable
reference observations remain available to audit earlier choices. This
migration requires the build's metadata and references, followed by bundling;
it cannot be performed on the immutable served bundle. Back up the build before
an explicit `--rebuild`, which intentionally discards historical observations
as documented by the authority CLI.

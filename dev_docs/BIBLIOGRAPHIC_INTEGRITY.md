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

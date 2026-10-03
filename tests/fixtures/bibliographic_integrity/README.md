# Bibliographic regression source values

`source.bib` contains selected fields from the library revision identified in
`provenance.json`; the fixture is about field preservation, not independent
validation of that library's catalogue. No PDFs or whole corpuscle are needed.
The source's `Dow. T.` author slot and article identifier stored in `pages` are
intentional: software must preserve supplied values without silently curating
them. The Delle Chiaje entries identify volumes in their titles and supply no
`volume` field, which the regression must not invent.

The recorded document hashes were read from retained bundle metadata. That
bundle is **v1.2.1**, built 2026-09-01 from pipeline SHA
`31784560ceaed157192ba2e4fa35a5f5c4b6feca`. It is older than the v1.4.0.dev0
bundle cited by the issues. Read-only inspection found:

- Mańko/Pugh metadata explicitly records BibTeX extraction and the entry key;
  the corresponding work retains correct Unicode authors but has a null
  `bib_imported_at`. This proves lost provenance in this retained bundle, not
  the later bundle's reconciliation history. This older authority schema has
  no reconciliation-decision table.
- Chun's document metadata retains the Excretionsporus title while its scalar
  authority row points to the Ctenophoren title. Current-code regression tests
  separately exercise rejection and membership rematerialization.

The tests materialize a fresh temporary authority from the pinned source
values and exercise formatting, graph roots, export/import and unchanged
refresh. `test_bibliographic_bundle_roundtrip.py` additionally runs the real
bundler and a separate production MCP server process, then calls
`format_citations` over stdio by document hash and work ID. It repeats that
boundary after explicit-import precedence and after export/import plus an
unchanged authority refresh. File hashes verify that formatting leaves the
bundle unchanged. Test output retains the bundle manifests, tool responses,
file hashes and server logs under pytest's temporary directory.

The prepared input boundary is parsed BibTeX metadata. There is no PDF
extraction, Grobid, embedding, taxonomy or figure acceptance in this test.
Church's deliberately conflicting imported locators are a labelled synthetic
precedence control; the real source entry is restored before source-value
assertions. Ahuja supplies its article identifier in `pages`, not a dedicated
`eid`/`articleno` field. Dedicated article-ID and chapter renderings have separate
synthetic coverage in `test_bibliographic_integrity.py`. The supplied Siebert
page range is preserved, not independently adjudicated or silently corrected.

The fresh Mańko bundle assertion establishes build-with-BibTeX provenance and
ordered Unicode authors. These tests do not mutate the retained bundle,
reconstruct the audited deployment's reconciliation history or establish that
the deployed v1.4 bundle has been rebuilt. Those remain distinct acceptance
requirements.

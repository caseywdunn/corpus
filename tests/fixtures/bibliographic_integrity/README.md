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
refresh. They do not mutate the retained bundle or establish that the deployed
v1.4 bundle has been rebuilt. That remains release acceptance work.

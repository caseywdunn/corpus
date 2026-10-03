# Hosia completed-gold figure acceptance, 2026-09-18

[receipt.json](receipt.json) records a measured one-document production replay:
completed gold artifacts → authority and taxon materialization → immutable
bundle → figure tools, live MCP and HTTP. It is an acceptance record, not a
replacement source fixture or a declaration that the entire candidate passed.
No figure evidence, captions, rights or derived database rows were manually
repaired for the replay. No extraction, OCR, model or embedding work ran.

The original PDF is Hosia et al. (2024), DOI
[10.3389/fmars.2024.1421514](https://doi.org/10.3389/fmars.2024.1421514),
source SHA-256
`1648cd91e97350ebd146034f9ab00da07ebe6c0617d0fdeb822f7fdcf544a413`.
The document completed successfully at `2026-09-18T19:30:01Z` under frozen
extraction commit `6fbf4e0`. Post-build, packaging and serving used `9d1eac1`.
Figure-rights, caption-taxon and figure-tool code did not change between those
revisions. The receipt preserves both exact commits and producer fingerprints;
its bundle manifest records the packaging revision, not the extraction revision.
Later surname v3 and other producer changes still require a refreshed candidate.

## What passed

| Issue | Measured result | Bound of the claim |
| --- | --- | --- |
| #302 | Figure 1 retains its explicit publication-license exclusion and page/caption provenance through the real bundle. Strict whole-image, labelled fallback, ROI, URL, direct HTTP and signed HTTP paths refuse; report permits delivery. Missing authority remains distinguishable from the recorded exclusion. | Hosia Figure 1 has no located ROI. Separate focused tests exercise actual excluded pixel crops and mixed-image inheritance. Figure 10 is permitted using the actual library license; this is software-policy evidence, not independent legal vetting. |
| #322 | Figures 4 and 5 exactly match the joins of retained source-checked caption fragments, with scale/anatomy text intact. Figure 1 retains the exclusion from its recovered fragments. | All retain `caption_completeness=unverified`. This does not establish universal completeness or the separate Sutherland case. |
| #323 | Exact caption span `[174,184]`, `N. septata`, resolves to taxon `1775295`, Nanomia septata, using caption-local Nanomia context. Caption-only discovery includes Figure 1 and retains Figures 10–13. Direct lookup, existing discovery tool and dossier return identical stored evidence; the materialized DB is copied byte-for-byte into the bundle. | Competing genera, punctuation, aliases and boundary negatives are separate focused tests. This does not estimate missed-abbreviation prevalence or grant figure rights. |
| #324 | Figure 12 preserves expected A–E, `no_labels_found`, no invented ROIs, and E returns whole-image fallback with its caption description. Figure 10 distinguishes located A from unlocated B. | Does not prove the named Siebert A–N/A–U or every panel description case. |
| #332 | Stored total, array length, manifest, bundle_info, corpus_summary and list_papers all report 15 records. | Agreement counts records, not unique scientifically correct figures, and is not a corpus-wide inconsistency census. |

The direct matrix has **27 cases**: three labels × three profiles for explicit
exclusion, inherited-permitted and deliberately absent-authority controls.
Its HTTP requests use in-process ASGI. Refused strict URL issuance is followed
by a locally minted valid capability to verify that the HTTP boundary independently
refuses delivery. The absent-authority control removes only the loaded database
handle; it does not alter source artifacts or rights.

The separate production stdio process passed **25 MCP calls and three actual
signed loopback HTTP downloads**. Each report download has the same hash as its
inline image. Inline strict refusals have MCP `isError=true`; JSON-returning
URL/ROI tools preserve their existing structured `forbidden` result contract.
The receipt records those distinctions rather than treating every refusal as
the same wire envelope. Both original document and bundle inventories stayed
unchanged; panel crop cache files were outside the bundle.

## Reproduction inputs and procedure

This compact record contains hashes and observed results, not the PDF, figure
images, complete artifacts, raw tool transcripts or executable test doubles.
The original raw replay-file hashes identify the measured run but cannot
reconstruct those files. Repeating this acceptance requires an actual completed
source artifact directory and taxonomy database matching the receipt's hashes,
or a separately built candidate whose changed identities are reported honestly.
Do not relabel a fresh rebuild as the historical frozen artifact.

Required inputs are:

- The original PDF identified above (worked-example library-relative name
  `H_J/Hosiaetal2024.pdf`) and completed `documents/1648cd91e973` artifacts.
  `source_document_artifact_sha256` pins every copied file, including images,
  original metadata, figures, chunks, summary and pipeline state.
- The actual taxonomy database, pinned by `packaged_artifact_sha256` and the
  caption-taxon input fingerprint. The original metadata includes the library's
  CC-BY-4.0 record; do not invent or override a license.
- The post-build/serving checkout and the existing
  [Hosia caption geometry fixture](../hosia_captions.json) for the independent
  retained-fragment comparisons. No model-dependent query is needed.

Reproduction sequence:

1. Verify PDF and completed-artifact identities, successful summary, and the
   extraction/materialization/cross-reference receipts. Copy the document and
   taxonomy into a new scratch build; record its byte inventory.
2. Call production `bib.authority.create_schema` and `phase1_corpus_papers` on
   its authority database. Call `pipeline.taxon_mentions.create_schema` and
   `build` on its mention database. Verify zero errors and unchanged copied
   document bytes. These calls generate real authority and caption links;
   do not manually insert derived records.
3. Call `mcpsrv.bundle.package` with PDFs omitted and no embeddings. Verify
   authority, taxonomy and mention DB bytes survive packaging, then open the
   bundle using `CorpusIndex`, `BiblioAuthority`, `TaxonomyDB`, `TaxonMentionDB`.
4. Replay `boundary_matrix.cases` using the actual figure tools,
   `make_figure_app`, bearer wrapper and production figure signer. Compare the
   compact recorded fields; ephemeral URL signatures are intentionally absent.
   For missing-authority cases temporarily omit the loaded authority handle.
5. Start `python -m mcpsrv.main <bundle> --transport stdio`. Send
   `live_transport.tool_calls[*].arguments` to their recorded tool names.
   Fetch the three report URLs through that process's loopback HTTP companion;
   compare downloaded bytes with inline image hashes. Hash the bundle again.

New packaging timestamps/database identifiers may differ in a new replay;
record new artifact hashes instead of rewriting this dated evidence. Compare
semantic contracts and original source identities. Stage fingerprint hashes in
the receipt use the exact recorded `json.dumps(..., sort_keys=True)` encoding;
the shared caption-evidence digest uses UTF-8, sorted keys, no ASCII escaping,
and compact JSON separators.

Focused regression controls reside in `tests/test_figure_rights.py`,
`tests/test_figure_licensing_states.py`, `tests/test_signed_figure_urls.py`, and
`tests/test_caption_taxa.py`. The durable receipt supplements those runnable
controls; asserting its saved booleans would not rerun production acceptance.

## Remaining observations

Hosia still reports Figure 2 missing, classifies the Table 1 image as a figure,
and contains an apparent cross-publication `24` in its missing-figure queue.
Source caption spacing artifacts remain. These facts prevent claims of complete
figure recall, correct classification or an adjudicated missing queue. The
one-document result is not complete gold/deployment acceptance, historical
v1.4 audit reconstruction, scientific-text recall, legal review, or retrieval
quality measurement. The stored source and replay identities remain fixed even
when later candidate builds improve those independent issues.

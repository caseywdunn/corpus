# Caption taxon evidence

The post-pipeline taxon index builds caption links from the final `figures.json`,
document `chunks.json`, and current taxonomy snapshot. This happens after the
separately scheduled figure passes, so changed captions do not inherit stale
annotation results. Per-paper receipts hash all three inputs and stamp the
caption policy. Unchanged builds skip writes; changed or removed figures,
document context, taxonomy and policy re-derive links.

Explicit full names and aliases resolve through the taxonomy with word
boundaries. Abbreviations reuse the taxonomic expansion rule: the genus must
be written in full in the caption or document, the expanded binomial must
exist in the snapshot, and surviving candidates must agree on one accepted
taxon. Caption context takes precedence over document context. An epithet by
itself never establishes a link. Competing genera and missing context remain
unresolved, with candidate names and a reason.

`taxon_mentions.sqlite` stores `caption_taxon_evidence` and
`caption_taxon_links`. Evidence includes the printed span, accepted taxon,
method, contextual genera, supporting document spans when applicable, and
input hashes. Responses cap displayed matches at 200 and unresolved entries at
50 per caption, with explicit totals and truncation flags. The figure lookup,
dossier, and direct figure record expose the same `caption_taxa` evidence.
Caption-only links can discover a paper even when its chunk taxon index lacks
the name. These taxon links do not establish caption ownership or image rights.

Legacy bundles remain readable with full-name-only caption filtering and
`availability: legacy_unavailable`; the server never expands abbreviations.
Rebuild the taxon index and bundle to enable contextual caption links. A
malformed update or a missing previously indexed figure artifact blocks build
success while retaining prior rows for diagnosis.

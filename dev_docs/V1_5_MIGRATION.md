# Updating a corpuscle for v1.5

This is the candidate upgrade procedure. Release acceptance, including the
full-corpus replay, remains pending in [PLAN.md](PLAN.md). The examples describe
the existing CLI; they do not declare an unverified migration complete.

v1.5 repairs stored bibliographic identities, extracted text and context, figure
geometry/captions/clearance, and the indexes built from that evidence. Install
the candidate on the build machine and rebuild from the library inputs before
replacing the served bundle. A server-only update can expose the new response
contracts but cannot supply evidence absent from an old bundle.

## Keep the inputs and comparison evidence

Retain the original PDFs, the configured `.bib`, lexicon and taxonomy source,
`config.yaml`, and any separately maintained curated BibTeX imports. Copy
curation changes back into the library so a clean build can reproduce them;
do not overwrite the library with an unreviewed export from a damaged authority
DB. Keep the previous build and served bundle as comparison/rollback artifacts.

Normal authority updates retain raw reference observations and decision history.
An explicit authority `--rebuild` discards that history. The top-level
`corpus run --force-rebuild` rebuilds cross-paper databases; it is not a way to
force every PDF through extraction. Neither flag is required for the normal
receipt-driven update below.

## Rebuild on the build machine

After installing the candidate with the normal [installation procedure](../INSTALL.md),
run from the corpuscle directory:

```bash
corpus --config config.yaml check
corpus --config config.yaml run --dry-run
corpus --config config.yaml run --no-bundle
corpus --config config.yaml status --report
```

Check that the configured OCR packs, Grobid service, embedding producer and
requested vision backend are available. Inspect the selected backend in the
run output: the OCR panel floor does not validate a configured vision build.
If the build runs in separate CPU/GPU phases, use the established phase
workflow and complete vision, embedding and post-processing before bundling.
On a GPU allocation, use `compute.accelerator: auto` or `require` with
`corpus run --require-gpu` so an unusable accelerator stops the run. An explicit
device setting is honored verbatim and bypasses that capability check.

The updated extraction policies and local producer identities invalidate
dependent work even when a PDF's hash has not changed. The normal run then
regenerates chunks, annotations, embeddings and cross-paper evidence as needed.
Allow for re-extraction, OCR and model inference; this is not just an SQLite
schema change. A model or curated surname-catalog change can also invalidate
consumers. See [surname recovery](SURNAME_RECOVERY.md) for its scope and costs.

For an independent clean candidate, copy the configuration and set its
`output_dir` to a new empty directory. Keep the source paths and build settings
equivalent. Run incremental checks in place: relocating an existing build tree
changes its embedded paths and is not a valid unchanged-rerun comparison.

## Check the candidate before bundling

Review structured stage failures and quality flags. A completed process does
not establish source fidelity, and an unresolved recovery receipt is not a
confirmed correction. Apply the project's [testing workflow](TESTING.md), plus
the reported source cases relevant to the corpus. Compare logical document,
figure, chunk and reference mappings with the retained build, and explain
changes rather than accepting matching totals alone.

Text changes require new embedding rows. Reference repair requires authority
materialization and reconciliation. Caption taxon links and figure-specific
clearance require their build-time evidence. A partial phase run cannot stand
in for these dependencies. For semantic retrieval acceptance, draw independent
source-review cases from the established fixture, then run the frozen queries
against the full reference and candidate corpora as described in
[retrieval evaluation](RETRIEVAL_EVALUATION.md).

Once the build and its acceptance checks pass:

```bash
corpus --config config.yaml run --only bundle
```

Serve that candidate bundle with the updated server and replay the relevant
MCP calls, including strict image/URL delivery. Validate its manifest and actual
stored counts. Replace the production bundle only after that review, following
the normal [deployment procedure](../DEPLOY.md).

## Client-visible changes

- A repaired bibliographic identity can change a `work_id`, including splitting
  previously conflated publication parts. Resolve saved document hashes again;
  do not assume an old work identifier still names the same publication.
- Re-extraction and chunking can change chunk and figure identifiers. Refetch
  current references before drilling into a saved result.
- Ordinary queries continue to read older bundles. New treatment/diagnosis
  filters return `rebuild_required` when their materialized context is missing;
  unavailable caption-link provenance likewise remains explicit.
- New scope, pagination, provenance and refusal fields explain what was
  selected, omitted or withheld. Source-text previews are bounded and can be
  truncated even when the requested chunk's own text is returned in full.
- Source-supported reference corrections affect derived mappings. The original
  parsed year, author and raw citation remain available for review. Ambiguous
  identities and unresolved OCR readings are not silently promoted to facts.

See the [MCP tool reference](MCP_TOOLS.md),
[bibliographic identity](BIBLIOGRAPHIC_INTEGRITY.md), and
[caption taxon evidence](CAPTION_TAXON_LINKS.md) for the individual contracts.

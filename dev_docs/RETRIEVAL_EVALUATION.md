# Source-graded retrieval evaluation

`tools/qc/retrieval.py` measures frozen query workflows against source-reviewed
passages. It calls the existing `get_chunks_for_topic` implementation and does
not change ranking, embeddings or artifacts. This is operator tooling; nothing
in `pipeline/` or the server imports it.

The [siphonophore worked example](examples/siphonophore_retrieval_evaluation.json)
preserves the seven exact requests and three prose controls from issue #320,
plus two explicit historical-intent controls. Diagnostic/key source targets
come from small fragments already used for extraction regressions. The three
prose controls now also have labels verified against original PDFs and rendered
paragraph regions in the [prose source review](examples/siphonophore_prose_source_review.json).
The base manifest has no sufficient independent sample. The separately frozen
[September 22 source-first supplement](examples/siphonophore_retrieval_review_2026_09_22/README.md)
adds ten source-reviewed questions across five papers while preserving all fixed
queries and the earlier insufficient materialized-unit draw. Its selection is
purposive and rank-unseen, not representative. Neither manifest has demonstrated
release retrieval acceptance; thresholds are predefined targets, not measured
quality. Use the supplement's exact frozen manifest for both full-population
captures.

## Freeze the evidence before tuning

Use the existing gold corpuscle as the **independent source-sampling
population**. A separate large fixture corpuscle is unnecessary. Generate an
independent set from its materialized diagnosis/key units before inspecting
their retrieval rankings:

```bash
python tools/qc/retrieval.py sample \
  --manifest dev_docs/examples/siphonophore_retrieval_evaluation.json \
  --output-dir /path/to/rebuilt/gold/output --count 20 --seed 3202026 \
  --out /tmp/retrieval-review.json
```

The sampler excludes the fixed manifest's target pages, deduplicates split
chunks describing the same source unit, and uses a seeded paper order so a
large paper cannot occupy the whole sample. It records the full selection,
eligible population count and digest, seed, and selected-paper count. This is
a sample of **materialized** diagnostic/key units in that corpus, not all
possible scientific questions or unrecognized treatments. Missing or small
populations remain explicit; do not silently replace difficult selections.

The sampling population is distinct from the **retrieval population**. Keep the
gold build identity with the sample's population digest and source review.
Run the frozen queries against the full intended served reference and candidate
corpora, including fixed audit papers outside the gold set. This preserves the
corpus-wide competition that produced table crowding and keeps the original
paper-scoped requests valid. Sampling from gold does not restrict retrieval to
gold or change any query's paper filter.

Review the selected PDF pages independently of retrieval results, write a
natural query about each source unit, and add source targets. The generated
queries are drafts; finalize them before capture. Source-grade any unreviewed controls as well. Each target records:

- A stable target ID, PDF short hash and physical pages.
- `source` evidence, including the source PDF hash and page/region or committed
  source-fragment pointer, with the review basis.
- `review_status="source_verified"` only after source review.
- `all_text` passage anchors and grade 0 (irrelevant), 1 (topical context), or
  2 (answer-bearing diagnostic/key/prose evidence).

Use distinctive content anchors, not a taxon name or heading alone. A key
heading in a table of contents does not establish retrieval of key branches.
The Erenna nectophore request already retrieved part of that key in the audit;
the manifest credits a relevant branch and separately tests the exact bract-key
request. A grade-2 hit means useful evidence was retrieved, not that a full
species comparison or complete identification key was recovered. Known-target
coverage is reported separately. Unmatched results remain **unjudged**, not
proven irrelevant; this is not precision or exhaustive recall measurement.

Freeze and retain the reviewed manifest before changing a ranker. Its complete
SHA-256 is recorded in each capture. The scorer refuses a capture from a
different manifest, preventing after-the-fact label changes from silently
changing the result. Correcting a label requires a documented new manifest and
both reference and candidate captures under that manifest.

## Source review of the worked-example prose controls

The source review on 2026-09-18 used retained `1.2.1` chunks only to locate
candidate paragraphs. Original library PDFs were then hashed, read and visually
checked as rendered page regions. No semantic retrieval ranks were inspected.
The companion JSON stores compact native excerpts, PDF hashes, physical pages,
boxes in PDF points, rendering settings and crop hashes. PDFs and raster crops
are not added to the fixture corpus.

| Control | Verified source passages | Grading scope |
| --- | --- | --- |
| Pneumatophore structure and gas gland function | Pugh 1983, physical p. 9; Wittenberg 1960, physical p. 6 | Gland location/ectodermal structure and carbon-monoxide production. Wittenberg separately qualifies the inference about elaboration of all float gases. |
| Erenna lures and fish prey | Haddock et al. 2005, physical p. 1; Haddock et al. 2017, physical p. 8 | Source observations and interpretation connect the appendages to fish attraction. General introductory luminescence prose is grade 1, insufficient alone. |
| Somatocyst function | Grossmann et al. 2014, physical p. 16; Pugh et al. 2018, physical pp. 8–9 | Lipid storage is discussed directly; proposed buoyancy/trim roles remain hypotheses. The reduced somatocyst in several diphyids is a source counterexample to a universal buoyancy claim, and relevant evidence for this query. |

These are known answer-bearing passages, not exhaustive relevance judgments.
The source authors' qualifications are part of the evidence. A topical mention
or an unjudged alternative is not silently promoted to a positive hit.
The nine source labels include eight grade-2 passages and one grade-1 context
passage; these counts describe labels, not observed top-k retrieval outcomes.

## Capture the unchanged and candidate workflows

```bash
python tools/qc/retrieval.py capture --manifest /tmp/retrieval-reviewed.json \
  --output-dir /path/to/full/reference/output --label retained-reference --role reference \
  --out /tmp/retrieval-reference.json
python tools/qc/retrieval.py capture --manifest /tmp/retrieval-reviewed.json \
  --output-dir /path/to/full/candidate/output --label candidate-build --role candidate \
  --out /tmp/retrieval-candidate.json
```

Capture runs the original call with its exact arguments and separate calls at
5 and 10 results. Corpus-wide requests remain corpus-wide; original paper
restrictions remain intact. The local route uses the bundle's query embedder
and the same MCP retrieval function, and may load its model. Use the normal
offline model cache for controlled runs. Retain raw results, bundle manifest,
embedding identity, paper-population digest and run label. Query failures are
recorded as operational errors and block acceptance.

Record the reference and candidate retrieval populations independently of the
gold sampling population. Check their paper inventories and explain any
membership difference before attributing a ranking change to the implementation;
the paper count alone is insufficient. All fixed audit papers and sampled source
papers must be present in both intended retrieval populations. A gold-only
capture can diagnose local behavior, but its hit rates do not establish
deployment retrieval quality or satisfy full-corpus release acceptance: removing
competing documents changes the task, even with identical query strings.

Do not identify a retained older output as the audited deployment. In the
worked example, `output1.2.1` is a separate reference; the issue audited
`1.4.0.dev0`, bundle timestamp `2026-09-09T09:23:59Z`, pipeline SHA
`734be4cfc301531bdcf90a1d23c2956bcbeed2fd`. Neither identity substitutes for a
fresh candidate build. Capture records `historical_audit_equivalence` as
unasserted; compare recorded identities explicitly.

To record assisted recovery, append a run with `mode="assisted"`, the same
`query_id`, a unique `variant` describing the workflow, its actual `call` and
returned `rows`. Paper discovery, narrowed queries and direct diagnosis fetches
belong here. These rows are scored separately and never improve an unassisted
gate. Preserve the actual assistance steps alongside the capture.

## Score and compare

```bash
python tools/qc/retrieval.py score --manifest /tmp/retrieval-reviewed.json \
  --capture /tmp/retrieval-reference.json --out /tmp/reference-score.json
python tools/qc/retrieval.py score --manifest /tmp/retrieval-reviewed.json \
  --capture /tmp/retrieval-candidate.json --out /tmp/candidate-score.json
python tools/qc/retrieval.py compare --reference /tmp/reference-score.json \
  --candidate /tmp/candidate-score.json --out /tmp/retrieval-comparison.json
```

`score` and `compare` exit 2 for failed or incomplete acceptance, while still
writing the report. The worked example's pre-tuning targets are:

| Group | Hit@5 | Hit@10 | Minimum source-reviewed set |
| --- | --- | --- | --- |
| Selected audit queries | 6/7 | 7/7 | All seven exact queries |
| Prose controls | 3/3 | 3/3 | All three controls |
| Historical-intent controls | 2/2 | 2/2 | Both controls |
| Independent diagnostic/key sample | 0.80 | 0.90 | Ten queries from at least five papers |

Missing labels, absent required runs, insufficient independent coverage, or
query execution errors block acceptance. Passing these gates is separate from
demonstrating an improvement: comparison requires increased independent hit
rate at either depth without a decrease at the other depth or in controls, and
a candidate that passes all gates. Historical evidence must remain retrievable
for historical requests; no universal recency preference is implied.

Each query/depth also reports unique-document diversity, the dominant-document
fraction, identified table rows, repeated rows from the same source table,
exact repeated table text, and known-positive-target coverage. Table
continuations may contain different relevant branches, so repetition alone is
not an error. Missing legacy table provenance yields an unknown repetition
rate plus the observed lower bound; it never becomes an assumed zero.
Source pages are checked when available. Old artifacts can match distinctive
source-reviewed anchors without page metadata, with that route counted
explicitly. Review unexpected misses and unjudged alternatives before making
scientific quality claims.

The selected audit queries are deliberately difficult examples. Their hit
rates and the sampled materialized-unit results must not be presented as a
deployment-wide failure percentage.

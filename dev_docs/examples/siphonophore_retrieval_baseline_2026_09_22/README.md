# Retained v1.2.1 fixed-query diagnostic baseline — 2026-09-22

The original read-only capture used the complete retained **v1.2.1** bundle and
MCP query implementation `d0014a34d48e45ff3fedb6f2f77ce09eec8e92ad`.
A normal guarded repeat at `8e11fee63dee06fef2f0452668d89a0836757863` now
records actual indexed-paper membership and table version 1775 before and after
all 66 queries. Its population guard passes. All 36 fixed calls returned identical
rows and distances, with identical metric rows and gates; no drift was found.
`guarded-receipt.json` and `guarded-command.json` identify this final baseline.
The original raw capture and score remain unchanged.
It is not the audited v1.4 deployment and is not a candidate/reference improvement
claim. The full 22-query manifest produced 66 serial calls; only the 12 fixed
queries and their 36 calls have been inspected. All independent ranks, distances,
metrics, aggregate gates and the full-score overall status remain withheld until
a candidate policy is fixed. No production ranking or labels were changed.

| Frozen group | Known-positive hit@5 | Known-positive hit@10 |
| --- | ---: | ---: |
| Audit | 1/7 | 2/7 |
| Prose controls | 2/3 | 2/3 |
| Historical controls | 1/2 | 2/2 |

These count source-verified frozen positive anchors. Unmatched passages are
**unjudged**, not proved irrelevant. In particular, several somatocyst results
appear to discuss the requested function but are outside the frozen passage
labels; their presence does not change this baseline's grades.

Observed fixed-query mechanisms:

- Nanomia's natural query returns Tung abundance-table text in 3/5 positions;
  the diagnostic formulation returns Tung in 6/10, and the German query in 4/5.
  Readable Latin names survive amid damaged text and repeated measurements.
- The paper-filtered exact Erenna bract-key request ranks its known positive
  chunks 155 and 156 at 7 and10. Its first five include a contents page,
  nectophore-key text, and other bract material. Document diversification alone
  cannot address this paper-filtered case. Erenna's broader audit question
  already hits a known nectophore-key passage at rank 2.
- The paper-filtered Apolemia lanosa request returns title/materials, two tables
  and two captions in its first five, while its indexed diagnosis chunk 25 is
  outside the first 10.
- The Physalia diagnosis target is indexed in Church chunk 53 but absent from
  its top10; returned historical material remains unjudged for that query.
  The historical controls must remain protected in any experiment.
- Table provenance is absent in this legacy bundle. Automated repeated-table
  rates are unknown (`null`), not zero. The abundance-table observations above
  come from inspecting the fixed returned text, not invented table identities.

A separate exact-availability scan used the evaluator's existing `target_match`
over **all complete retained chunks** of fixed target papers and then checked
actual indexed text. All **4,266 rows in 12 fixed target papers** are indexed,
with zero missing IDs and zero text differences. Every fixed query has an
available indexed positive. One positive Erenna cornuta bract branch, shared
by audit_05/audit_07, fails the single-chunk target match because its character
anchor occurs in chunk 154 and endpoint in 155; the remaining four bract targets
are available. No absent indexed positive explains any fixed query's hit miss.
This is exact anchor availability, not exhaustive semantic fidelity grading.

The original guard checked **artifact membership**: 1,775 papers and all 17
source-target papers were present. Artifact population digest:
`21e547692636fef6cca1cbfab9f44104691004fb9403078d950908e47ce03a0f`.
A subsequent projected full-index census found **1,774 indexed papers** with
259,952 rows, digest
`bdebea69fa74346e2afa29b43817fefbfdaf3b3136f852dd8c03b668076953c4`.
The sole artifact-only paper, Wangersky_Lane1960 (`62e07c061591`), has zero
materialized chunks and an empty-text quality error. There are no orphan index
papers. Thus artifact membership and actual retrieval competition must remain
separate; the original artifact-only guard was insufficient.

`indexed-population-census.json` preserves the full per-paper row counts,
reconstructible sorted inventory and evidence for the empty artifact. Its version
checks (1775 before and after) bracket this later census, **not** the original
queries. It is diagnostic evidence; it does not retroactively add query-bracketing
checks to the original capture. The later **normal** guarded capture separately
records version 1775 around all queries and reproduces the exact same full census.
Its compact receipt references that identical stored count map instead of
duplicating 1,774 entries. All 17 target papers are indexed; no paper population
was removed or added between the two baseline runs.
The retained document embeddings expose BGE-M3/1024 dimensions but no versioned
producer receipt. This limitation precludes presenting the run as a producer-
verified embedding migration or splicing in new document embeddings.

Original capture cost: 43.54 seconds wall; peak RSS 3,663,600 KiB (about 3.49 GiB);
zero swaps. The guarded repeat cost 47.06 seconds, peak RSS 3,647,528 KiB (about
3.48 GiB), zero swaps and aggregate CPU 301%.
The requested offline/CPU/OMP/MKL/OpenBLAS/Rayon settings are preserved in
`command.json`. Original native aggregate CPU was 305%, despite thread variables set to 1; this
was not literally a one-core process. No OCR, document embedding, network request
or source-library mutation was performed. The model process has exited.

Reproduction and durable evidence:

- `command.json` / `receipt.json` preserve the original run.
  `guarded-command.json` / `guarded-receipt.json` preserve the final normal
  guarded run, exact implementation/environment, raw output hashes, actual
  index membership and query-bracketing versions, fixed-only drift checks,
  and measured resources. Raw captures/scores contain withheld independent
  outcomes and are intentionally not copied into this directory.
- `fixed-results.json` contains every fixed query's original/at5/at10 metric row
  and the six fixed group gates. It contains no independent metrics or rows.
- `target-availability.json` retains source anchors, matching chunk IDs, complete
  artifact hashes and matched passage hashes without repeating passage prose.
- `index-consistency.json` retains the read-only scalar filter, table version,
  complete-paper row/text checks, all fixed target matches and the selected
  index-row digest. It uses no embedding or vector search.
- `guarded-validation.json` records reconstruction, checksum and fixed-only
  comparison checks. `SHA256SUMS` protects these durable projections. Old absolute scratch paths
  are provenance, not portable dependencies. The retained bundle/source build
  is external and must remain available to repeat the measurement.

To reproduce the complete baseline, check out the exact implementation in
`guarded-command.json`, provide that same immutable v1.2.1 bundle and use the frozen
manifest in the sibling `siphonophore_retrieval_review_2026_09_22` directory:

```bash
env HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 CORPUS_DEVICE=cpu \
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 \
  TOKENIZERS_PARALLELISM=false python tools/qc/retrieval.py capture \
  --manifest dev_docs/examples/siphonophore_retrieval_review_2026_09_22/manifest.frozen.json \
  --output-dir /path/to/retained-v1.2.1/_serve \
  --label retained-v1.2.1-full-guarded-diagnostic-20260922 --role reference \
  --out /tmp/retained-capture.json
python tools/qc/retrieval.py score \
  --manifest dev_docs/examples/siphonophore_retrieval_review_2026_09_22/manifest.frozen.json \
  --capture /tmp/retained-capture.json --out /tmp/retained-score.json \
  > /tmp/retained-score-command.json
```

The score file contains the full report; while independent outcomes are sealed,
inspect only a projection filtered by manifest `group != "independent"`. Do not
inspect the overall status or infer it from the score command exit code, because
it combines independent outcomes. Operational call-error counts may be checked
separately.

The availability scan selected only those 12 fixed queries, read all chunks of
their positive-target papers, attached `paper_hash`, and applied evaluator
`enrich_rows` then `target_match` without changing either function. The separate
index read used `table.search().where(predicate).select(["text", "metadata"])`
with the recorded scalar predicate and a limit one greater than the complete
artifact row count. It required no excess/duplicate rows, checked every chunk ID
and exact text against the artifacts, and repeated the same target matcher on
indexed text. This checks materialization and indexing separately from ranking.

No candidate result has been measured. Existing source-text/context corrections
should first be evaluated against an appropriate comparable candidate before
considering additional ranking behavior. A different paper population, a partial
thesis replacement or an unversioned embedding splice would not be a valid
same-background comparison.

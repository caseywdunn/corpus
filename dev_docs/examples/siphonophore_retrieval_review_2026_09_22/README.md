# Separately frozen source-first diagnostic review (#320)

This is a **purposive, retrieval-rank-unseen supplement**: ten distinct source
units from five papers, with one natural query per unit and sixteen known
answer-bearing passage anchors. It is not a probability sample or a deployment
quality result. No query embedding, semantic retrieval, OCR, extraction or
ranking experiment was run for this review.

The earlier [materialized-unit sample](../siphonophore_retrieval_review_2026_09_18/README.md)
remains unchanged and insufficient: five selected units/two papers, of which
three units/one paper were source-scorable. Its false-key findings and limited
coverage are not replaced by this supplement. The saved gold artifacts all
carry treatment policy v2; even read-only v3 treatment-context replay leaves
resolved diagnoses confined to Totton. Source selection here therefore does
not require a successful materialized diagnosis tag.

`protocol.frozen.json` fixes the five papers and ten candidate units before
source-image grading and query formulation. Paper selection followed saved
source-text inspection, **not random sampling**. Four Vanhöffen units are the
first four numbered species treatments in source order. Totton's two selected
species keys differ from the earlier reviewed genus keys. All candidates were
retained; none was replaced to reach ten. The selection excludes the fixed
manifest's source target pages and every page in the earlier frozen draw.

`source-review.json` records direct visual inspection of all seventeen original
PDF page renders, source hashes, precise region boxes, crop hashes, relevance
judgments and qualifications. Review was by Codex agent `/root/fix_figures`;
a separate Codex agent `/root/fix_bibliography` checked only Vanhöffen pp. 7–8.
This is agent review, not blind, human or domain-expert adjudication.
Historical names and tentative synonymies are retained as the sources state
them; the review makes no current-taxonomy claim.

| Unit | Paper; physical source pages | Diagnostic content |
| --- | --- | --- |
| 01 | Totton 1965a; 58 | Agalma elegans/haeckeli bract characters in a species key |
| 02 | Totton 1965a; 164 | Lensia subtilis/meteori somatocyst and pedicle characters in a species key |
| 03 | Totton 1965b; 3 | Lensia baryi anterior nectophore compared with L. achilles |
| 04 | Totton 1965b; 4–5 | Separate L. cordata treatment: somatocyst shape/position and ridge crests |
| 05 | Margulis 1984; 3–4 | L. campanella elongata subspecies diagnosis and original English summary |
| 06 | Gasca and Suárez 1993; 1 | L. canopusi contrasted with L. hotspur and L. cossack |
| 07 | Vanhöffen 1906; 4–5 | Sphaeronectes gracilis/Köllikeri reservoir orientation |
| 08 | Vanhöffen 1906; 5–6 | Muggiaea atlantica/kochi hydroecium and reservoir extent |
| 09 | Vanhöffen 1906; 7 | Galeolaria truncata nectophore morphology and reservoir extent |
| 10 | Vanhöffen 1906; 8 | Separate G. biloba treatment: hooked teeth on the bract margin |

The units are taxonomic treatments or keys, not extra questions about the same
fact. Two source passages from one treatment remain one query. A grade-2 target
can supply one useful branch or part of a comparison: a hit is not proof that
all diagnostic characters or both species were retrieved. Known-target coverage
is reported separately by the evaluator. Unmatched content remains unjudged.

`manifest.frozen.json` appends these ten queries to the twelve unchanged fixed
audit/control queries and their existing labels. It preserves all predeclared
hit@5/hit@10 thresholds. Freeze identities are in `receipt.json` and
`SHA256SUMS`; the evaluator uses a canonical JSON digest, distinct from the
file-byte SHA-256. The original sample's byte hash is pinned too.

## Reproduce source regions

PDFs, full Docling captures and raster images are deliberately not vendored.
Absolute `/tmp` and library paths in frozen JSON are historical provenance, not
portable dependencies. `protocol.frozen.json` supplies library-relative source
paths. With the same original PDFs and recorded PyMuPDF version, render the
full-page and smaller review crops into a separate directory:

```bash
python dev_docs/examples/siphonophore_retrieval_review_2026_09_22/render_review.py \
  --library /path/to/siphonophores/library --out /tmp/retrieval-source-review
```

The command checks all PDF hashes first and reports whether rendered image
hashes match. It renders existing sources only; it does not run OCR or modify
any library file. A different renderer version may change PNG bytes; report
that difference rather than changing the frozen expected hashes.

## Remaining acceptance

Source-label preparation meets ten queries/five papers for this separately
identified supplement. **Retrieval acceptance and improvement are unmeasured.**
Use [the existing evaluator](../../RETRIEVAL_EVALUATION.md) with this exact
manifest for both captures, after verifying the intended reference/candidate
paper inventories and model identities. All original and sampled source papers
must be present. The available older full bundle is a separate reference, not
the audited September 9 deployment. A gold-only run changes corpus competition
and cannot establish deployment retrieval acceptance.

The fixed audit, prose and historical controls, independent hit-rate and
improvement thresholds, repetitive-table/diversity measurements, and review of
unexpected misses/unjudged alternatives still apply. This selected supplement
cannot estimate failure prevalence, exhaustive relevance or retrieval recall
across the whole corpus. Later label corrections require a new manifest identity
and captures of both workflows; no ranking-guided edits are authorized here.

Validation: the existing evaluator accepts the manifest; all twelve fixed queries
and targets compare equal to the base manifest. Thirteen evaluator tests pass.
The portable render command reproduces all 32 expected hashes (17 full pages
and 15 smaller regions). See `validation.json`; these are artifact checks, not
retrieval outcomes.

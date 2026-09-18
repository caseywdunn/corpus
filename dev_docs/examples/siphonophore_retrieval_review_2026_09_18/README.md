# Frozen independent retrieval source review, 2026-09-18

The requested sample of 20 units (seed 3202026) yielded **five eligible units from two papers** in the completed 35-document gold build. Original-PDF review found **three scorable queries, all from one paper**, plus **two false identification keys**. Two usable units also carried incorrect materialized treatment names. The predefined minimum of **ten queries across five papers remains unmet**. No selections were replaced.

| Selection | Original source | Completed decision |
| --- | --- | --- |
| independent_01 | Ahuja et al. 2026, physical p. 25 | Bibliography entries 103–107, not an identification key. Unscorable; independently confirmed by a second Codex agent. |
| independent_02 | Totton 1965, physical p. 120, printed p. 115 | Prayinae genus key. Scorable radial-canal distinction; stored `Maresearsia praeclara` subject is incorrect here. |
| independent_03 | Totton 1965, physical p. 57, printed p. 52 | Agalmidae genus key. Scorable diagnostic branches. |
| independent_04 | Totton 1965, physical p. 153, printed p. 148 | Nematocyst length/diameter table, not an identification key. Unscorable; no organism inferred from nearby captions. |
| independent_05 | Totton 1965, physical p. 130, printed p. 125 | Prayoides diagnosis, confirmed by original genus heading on p. 129. Stored `Rosacea cymbiformis` subject is incorrect here. |

`sample.frozen.json` is a byte-identical copy of the original sampler output, including its fixed queries, five draft queries and exact selection. `reviewed_decisions.json` adds the previously completed review decisions, source-grounded natural questions and compact positive anchors for the three usable units. It does not replace or edit the original sample and is not itself a capture manifest. Both false keys remain explicit, with null natural questions and no positive targets. A later evaluation manifest must preserve these provenance and insufficiency facts rather than silently dropping them.

`receipt.json` pins the 6fbf4e0 extraction, a3e638d sampler, all 35 source PDF identities and snapshot artifact hashes, population/selection digests, review-decision hashes and source-freeze completion scope. Each review record retains the original selected source item/region and exact source/render hashes. Agent visual review is identified explicitly; this is not human or domain-expert adjudication. Source PDFs, page rasters, full extracted documents and absolute local paths are omitted from this compact handoff.

Original source pages were rendered with PyMuPDF 1.28.0 at 2×, without OCR. Native PDF text was used only to locate passages. Adjacent original p. 129 establishes the Prayoides genus context; it is not an additional sampled unit. Individual grade-2 anchors establish useful answer-bearing evidence, not complete key recovery or exhaustive relevance judgments.

The gold process wrapper returned exit 143 while its log reported all 35 documents processed and the extraction step successful. Before selection, the source-freeze operation separately verified successful document summaries and completed extraction/chunk/cross-reference receipts. This record retains the wrapper discrepancy and does not claim a clean process exit.

This consolidation only serializes already completed review work. No new source review, sampler modification, producer investigation, tests, retrieval calls, OCR or model work ran. False-key and wrong-treatment producer mechanisms remain uninvestigated. No hit rates or ranking improvements have been measured. The frozen 6fb gold population is distinct from both a later producer refresh and the full intended served reference/candidate corpora required for retrieval evaluation.

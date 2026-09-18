# Frozen selected spacing review (#334)

This compact record preserves the 13 source decisions frozen on 2026-09-18
before visual grading. It is follow-up acceptance evidence, **not issue or
release closure**. No corpus-wide precision, recall or detector-error total is
claimed. The reviewer was agent `fix_figures`; agent `fix_query_contracts`
independently reviewed only Bernstein's unresolved authority boundary. This was
agent visual review of selected known cases, not blinded or human labeling.
Prior preliminary Bernstein review is disclosed in the frozen population.

`population.json`, `grades.json`, `crops.json` and `independent-review.json`
are byte-for-byte copies of the original review records. Their historical
absolute paths identify the original captures and checkout; those paths are
provenance, not dependencies of the renderer below. The population and grades
were **not revised** after the later mixed-font proposal/warning fix.

| Stratum / row IDs | Whole-decision result | Boundary result |
| --- | --- | --- |
| Captured DuClos/Hissmann inputs replayed with retained evidence: P1, P2 | 2 correct / 2 | 15 correct / 15 |
| Actual Bernstein frozen-`6fb` gold repair: P4 | 0 correct, 1 unresolved / 1 | 1 correct, 1 unresolved / 2 |
| Available produced-repair records combined | **2 correct, 0 incorrect, 1 unresolved / 3** | **16 correct, 0 incorrect, 1 unresolved / 17** |
| Controlled reconstruction of reported Mapstone p68 collapse: P3 | 1 correct / 1 | 3 correct / 3 |
| Additional native DuClos geometric proposals: G2–G4 | 3 correct / 3 | 25 correct / 25 |
| Unchanged or rejected controls: N1–N5 | 5 correct / 5 | Not a repair-precision population |

The population includes all five long runs in the retained full DuClos p6
capture. G1, `werecollectedatFridayHarborLaboratoriesbetween`, has visible
source spaces but received no geometric proposal at the reviewed code snapshot.
It is recorded as an unresolved named recovery gap, not a recall measurement.
The subsequent mixed-font followup can propose its internal boundaries but
still leaves the phrase unchanged because regional OCR modes disagree; see the
[parent fixture notes](../README.md). That later evidence does not alter these
frozen labels.

G2–G4 have source-supported geometric boundaries, but **their individual
historical OCR decisions are unavailable**. The earlier aggregate observation
“two accepted, one disagreement” does not identify which proposal was accepted
or retain those OCR outputs. Neither this review nor this durable copy supplies
that missing evidence. These candidates therefore do not enter an historical
accepted-repair precision denominator. Later regional OCR is new evidence and
cannot retrospectively fill the missing historical decision receipts.

P4's genus/species boundary is visible; the additional space before `O.F.M.`
is too tight to judge confidently. The included [original-source crop](P4-source.png)
allows review without the original PDF. Both agents judged that printed
whitespace boundary unresolved, **not an established incorrect repair**. Do
not exclude P4 and describe the produced-repair population as 100% correct.

The five controls retain distinct limits: N1 shares P3's source region and is
already spaced; N2's closed Mapstone p200 phrase is actually printed without
spaces; N3 is the intact morphological word `nectophoral`; N4 is the genuine
19-letter German word `Einzelinformationen`, below the repair's 20-letter floor;
N5 is the rejected OCR-only `cou nt` split. They are not a random false-positive
population. Long-run screening hits are not automatically extraction errors.

## Re-render the source regions

The original PDFs are not vendored. From a checkout, using the Corpus Python
environment with **PyMuPDF 1.28.0** and the matching original library:

```bash
python tests/fixtures/table_structure/spacing_review_2026_09_18/render_crops.py \
  --library /path/to/siphonophores/library \
  --output /tmp/spacing-review-crops
```

Add `--rows P4` to render only the unresolved Bernstein region. The helper
verifies each full original PDF SHA-256, renders the saved physical-page,
top-left rectangle and DPI, then verifies the exact PNG hash. It never runs
OCR, a model, pipeline extraction, or grading. It repeats the original native
line inspection before rasterization to reproduce MuPDF's font resources.
Coordinates already include the original
crop padding. All 13 rows correspond to 12 PNG files because N5 reuses P1's
crop. Source hashes are in `population.json`; crop bounds and hashes are in
`crops.json`. All five source PDFs must match for the default complete render.

To verify this durable package without PDFs:

```bash
cd tests/fixtures/table_structure/spacing_review_2026_09_18
sha256sum -c SHA256SUMS
```

Original immutable review hashes:

- Population: `5e179cd3a822f317cc86e6a5e4c467557ab4b75f523a1641bc878b08e26d7c4b`
- Grades: `a646596d18ca774f2c2ebb01b6cc33ca96519701e9178d102b16abfdd80c135a`

Full Docling captures, full build receipts, source PDFs and duplicate source
rasters are intentionally not included. The frozen records preserve the
observed strings, proposed boundaries, grades, source identities and measured
limits, and permit independent rendering of every reviewed source region.
They do not reproduce missing model decisions or replace the retained
[source-fragment fixtures](../README.md).

Bernstein authority whitespace, G1 recovery, missing historical DuClos OCR
decision identities, named downstream annotation/embedding regeneration and
final-candidate whole-corpus acceptance remain separate outstanding gates.

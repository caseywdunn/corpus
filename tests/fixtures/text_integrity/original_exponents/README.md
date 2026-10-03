# Original scan exponent evidence (#303)

These compact observations come from `S/Stepanjants2014.pdf`, physical page 8.
`capture.json` pins the complete source PDF SHA-256, original source glyphs,
rendered crops, the OCR producer and outputs, and the actual prepared text item
from the retained `6fbf4e0` gold build. No PDF or built corpuscle is committed.
The prepared PDF and Docling artifact identities are recorded separately.

An agent visually reviewed the original rendered region before this repair.
Both printed expressions have a raised `3` after Cyrillic `мм`. This was an
agent source review, not a human/domain review. The hidden original text layer
encodes both digits as `3` with `raised=False`; normal full-page OCR subsequently
turned them into punctuation. The retained Docling item has `2000  мм?` and
`0.3-10 мм'`.

The source-only pilot ran on 2026-09-22 with Tesseract 5.5.2 and the recorded
English model. It scanned the 14-page original's character geometry and found
exactly two candidates, both on page 8. Only three bounded crop OCR calls ran:
PSM 7 and 10 independently read the first digit as `3`; PSM 7 returned empty
for the second, stopping that check. The PSM 7 disagreement is preserved and
the second expression remains unchanged. No expected-power lookup, whitelist,
dictionary, OCR prompt, or unit-specific digit substitution was supplied.

Both candidates also retain independent source-ink bounds for the adjacent
baseline letter and raised digit. Baseline placement is checked against actual
ink, because hidden OCR font geometry alone can be wrong. The first candidate
was then applied through production code to the **actual retained full gold
Docling document and prepared PDF**: item `#/texts/125`, character span
`[734, 735]`, changed from `?` to `³`; the second was recorded unresolved.
Its source and prepared PDF hashes are separate, and the raw `orig` survives.

The source recovery runs before PDF preparation can erase this evidence. The
normal extraction hook consumes its immutable receipt only if the prepared
page size, unique glyph location, numeric-unit prefix and structured local
text agree. A conflicting digit, missing/ambiguous owner or alignment refusal
does not modify text. A new nested producer identity invalidates preparation
and its downstream consumers, including when the digit OCR installation changes.

## Regression and limits

`tests/test_source_exponents.py` replays the exact retained crops and OCR
observations without invoking OCR. It also tests the actual parser/candidate
and source-ink logic, preparation receipt retention, immutable source evidence,
save/reload, materialized chunk provenance, idempotence, finite raster budgets,
and producer invalidation. Captured fragment tests use an isolated text item;
their local `#/texts/0` is not the full source artifact's `#/texts/125`.

Synthetic negative controls explicitly cover baseline/subscript/full-height
ink, bogus raised flags, ordinary affiliation/year text, a quantity on another
line, conflicting numeric values, distant prepared text and ambiguous owners.
These controls are not additional independently graded source examples.

No fresh full-page OCR, Docling model extraction, embedding or full candidate
bundle was run for this patch. Replaying verified original evidence onto the
retained prepared artifact is distinct from a fresh normal build. The source
receipt hook still needs that targeted refresh in the release acceptance run.
The candidate grammar covers encoded `2`/`3` in bounded Latin/Cyrillic numeric
unit contexts; it is not a general mathematical notation repair. The second
Russian exponent and the broader independently graded scientific precision/
recall acceptance remain open under #303.

To reproduce the source-only inspection with the original library available:

```python
import hashlib
import json
from pathlib import Path
import fitz
from pipeline.source_exponents import inspect_source_exponents

fixture = Path("tests/fixtures/text_integrity/original_exponents")
capture = json.loads((fixture / "capture.json").read_text())
source = Path("<library>") / capture["source_pdf"]
assert hashlib.sha256(source.read_bytes()).hexdigest() == capture["source_sha256"]
with fitz.open(source) as pdf:
    receipt = inspect_source_exponents(pdf)
print(json.dumps(receipt, ensure_ascii=False, indent=2))
```

Tesseract and its English model must be available on `PATH`; otherwise the
unchanged source candidates correctly remain unresolved. Producer differences
must be recorded when comparing a later reading with this frozen capture.

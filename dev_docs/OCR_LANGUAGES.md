# OCR language packs: what the gold set says

Everything here is measured against the 35-document gold transcription set
(761 pages, 1594–2026, 13 languages), whose pages were transcribed from
rendered images alone. Numbers are token coverage — the fraction of the words
actually printed on the page that came back.

Two measurement notes, because both cut against the conclusions:

- Absolute values in the per-pack tables come from `ocrmypdf` → `pdftotext`,
  not from the pipeline's docling path, so they are not comparable to
  `tools/qc/fidelity.py` output. Comparisons *within* a row are valid; across
  tables they are not.
- The scorer normalises away diacritics, so a model that renders `fran` for
  `från` is not penalised. That biases every result below **in favour of
  `eng`**, and `eng` still loses where the tables say it loses.

## 1. Availability dominates every other choice

19th-century German set in Fraktur, same pages:

| | coverage |
| --- | --- |
| poppler's embedded text layer | 0.05–0.15 |
| OCR with `deu_latf` | **0.80–0.83** |

Everything else in this document is worth 0.01–0.08. Pack *selection* is
tuning; pack *availability* is a cliff. `tools/install_tessdata.sh` exists for
this reason, and `scan_detection.json` records `tesseract_pack_available` so a
missing pack is visible rather than inferred from a bad score.

## 2. Models compete for the same glyphs

A second pack pays off when it covers something the first cannot, and costs
when it duplicates the first's territory. Both directions are large:

| document | packs | coverage |
| --- | --- | --- |
| `Kawamura1911a`, vertical pages | `jpn_vert` | **0.574** |
| | `jpn` | 0.246 |
| | `jpn_vert+jpn` | **0.186** |
| | `jpn_vert+eng` | 0.176 |
| `Linnaeus1735` (Latin throughout) | `lat+eng` | **0.624** |
| | `lat` | 0.562 |
| | `eng+deu+deu_latf+fra+lat+spa+por` | **0.552** |
| `Tilesius1814` (Swedish + Latin diagnoses) | `swe+cat+fra+eng` | **0.868** |
| | `swe` | 0.837 |
| | `swe+eng` | 0.833 |

`jpn_vert+jpn` is the extreme: one language, one script, two models contesting
identical glyphs — worse than either alone by 0.39. Seven packs on a
monolingual Latin text is the same failure, milder. Four packs on a genuinely
bilingual page is a gain.

`deu+deu_latf` is the constructive case and looks like a counterexample until
you see why: one language, two *letterform systems*, so the models do not
contest the same glyphs.

**The rule is not "more" or "fewer". It is: add a pack for content that is
actually on the page.**

## 3. Why the native pack matters: character coverage, not word lists

The obvious hypothesis is dictionaries. Tesseract exposes
`language_model_penalty_non_dict_word` (0.15) and loads six dawgs per language,
so a plausible story is that a mismatched word list mangles vocabulary it does
not contain.

**Measured, that story is wrong.** Re-running `Candeias1932` with
`load_system_dawg 0`, `load_freq_dawg 0`, `load_unambig_dawg 0`,
`load_bigram_dawg 0`:

| packs | dawgs | coverage |
| --- | --- | --- |
| `por` | on | 0.954 |
| `por` | **off** | 0.954 |
| `eng` | on | 0.825 |
| `eng` | **off** | 0.825 |

Not merely similar — **zero of 2,236 recognised tokens changed.** Those
parameters govern the legacy word permuter, and Tesseract 5 runs the LSTM
engine, so they are inert.

What actually differs is the **character repertoire the LSTM was trained on.**
Of the tokens `por` recovers and `eng` loses on the same pages:

| | tokens | carrying a diacritic in the original |
| --- | --- | --- |
| recovered by both | 964 | 110 (**11.4%**) |
| **lost by `eng`** | 158 | 132 (**83.5%**) |

A sevenfold enrichment. The losses are `dimensões`, `campânulas`, `contribuição`,
`região`, `nectóforo`, `identificação` — `eng` does not merely drop the accent
(the scorer would forgive that); it misreads the accented glyph as a different
character, so the token fails to match even after normalisation.

That predicts native-pack advantage should track diacritic density, and it does:

| document | diacritic letters | native | `eng` | native advantage |
| --- | --- | --- | --- | --- |
| `Linnaeus1735` (lat) | 0.3% | 0.562 | **0.600** | −0.038 |
| `DeHaan1827` (nld) | 0.4% | 0.789 | 0.785 | +0.004 |
| `Tilesius1814` (swe) | 2.9% | **0.837** | 0.751 | +0.086 |
| `Candeias1932` (por) | 3.4% | **0.931** | 0.850 | +0.081 |

Dutch is near-unaccented and `eng` matches it. Latin is unaccented *and* has a
weak traineddata, so `eng` beats it outright. Swedish and Portuguese are
diacritic-dense and the native pack wins by ~0.08.

## 4. Adding `eng` usually pays; substituting it does not

Because Tesseract arbitrates per word, a second Latin-script model beat **both**
its members in three of four documents: `por+eng` 0.944 over `por` 0.931 and
`eng` 0.850; `nld+eng` 0.818 over `nld` 0.789 and `eng` 0.785.

But the asymmetry between adding and replacing is severe. Splitting tokens by
English-wordlist membership — which here is a proxy for *language*, not for
dictionary effect — replacing `por` with `eng` costs:

- words in the English list: **−0.010**
- words outside it: **−0.129**

Thirteen times the damage, on the partition holding **89–92% of this
literature's tokens**. Historical systematics is binomials, anatomical terms and
non-English body text; that is where a mismatched pack does its damage, and it
is also the vocabulary retrieval depends on.

**Add a complementary model. Never substitute one.**

## 5. Pinning is not automatically better than detection

Pinning `ocrlang` on 31 of 35 gold documents from an annotated `doclang` moved
corpus-wide prose coverage **0.9474 → 0.9450** — 31 pages better, 53 worse. The
languages were right. Each pin replaced a detected *union* with a singleton.

Detection's failure mode is over-union (seven packs on a monolingual Latin
text). A pin's failure mode is narrowness. Neither is safe by default:
detection's union won on `Tilesius1814` (0.868 vs 0.833) and `DeHaan1827`
(0.823 vs 0.818); the pin won on `Linnaeus1735` (0.624 vs 0.552).

**Pin what detection cannot infer. Do not pin merely to be explicit.**

## 6. What a pin can express — and the one thing it cannot

| only a pin knows | detection handles it |
| --- | --- |
| Fraktur — `de-Latf` → `deu_latf`; langdetect can only ever say `de` | dominant body language |
| script variant — `zh-Hans` vs `zh-Hant` | script family |
| Ancient Greek — `grc`, which has no ISO 639-1 code | **writing direction** (geometric, #196) |

Writing direction is the trap. BCP-47 describes language and script; vertical
setting is typesetting, and is deliberately outside `bcp47_to_tesseract`'s
scope (#215) — so any pass deriving `ocrlang` from `doclang` emits `jpn` for a
vertically-set document **every time**. The pin then overrides the geometric
detector that had it right, and the hint that would have reported it was
suppressed precisely because a pin existed.

Found in the siphonophore library on 7 of 28 CJK documents, each worth −0.33 on
its vertical pages. `scan_detection.json` now records
`ocrlang_overrides_vertical_cjk` when a pin contradicts the measurement.

## Summary

Choosing OCR languages is not "identify the language". It is choosing a small
set of models whose character repertoires cover what is printed, without two of
them contesting the same glyphs — over a vocabulary that is mostly outside
every one of their dictionaries, which turns out not to matter.

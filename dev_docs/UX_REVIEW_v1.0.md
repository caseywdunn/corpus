# Pre-release UX review — clean install + smoke-test build

**Branch:** `ux-smoketest-review` (from `dev` @ 7a20d80, `corpus 1.0.0rc1`)
**Date:** 2026-08-01
**Method:** new-user persona working from `README.md` alone (no code reading to decide *what*
to do), on `~/data/siphonophores_smoke` — 32 PDFs, 676 pages, 13 languages, 4 Fraktur, 5
Russian, 3 image-only, one 314-page monograph.
**Host:** linux-x86_64, 12 cores, 31 GB RAM, no usable GPU (NVML driver/library mismatch),
Grobid 0.8.1 already running on :8070.

## Outcome

The build succeeded end to end: **55 minutes, exit 0, 32 papers / 3,272 chunks / 360 figures**,
every stage at 100% in `corpus status --report`, all four cross-paper artifacts present, and
`corpus serve --check` green. An MCP client answered substantive questions correctly, including
synonymy-resolved taxon queries and formatted citations.

The install and the CLI are in good shape. **The serious problems are all silent quality
degradation that the tooling reports as success** — three separate mechanisms where a run
exits 0, `corpus status` is all green, and the corpus is materially worse than the docs
promise. Those are P1 below.

---

## P1 — silent quality loss reported as success

### 1.1 OCR is skipped on 29 of 32 papers, including all four Fraktur ones — **FIXED**
The headline multilingual-OCR capability never runs on the material it exists for.

> **Resolution (2026-08-02).** CWD confirmed the intent: OCR everything that is not digitally
> native, even if it has been OCR'd before, because third-party OCR quality is highly variable.
> That was not happening and could not have been — `detect_scan_type` only ever inspected the
> *content* of the text layer, so `born_digital` meant "reads plausibly", not "produced
> digitally". See §1.1a below for what shipped.

`scan_detection.json` across the corpuscle: **29/32 `born_digital` / `clean_text_layer` /
`needs_ocr: false`**; 3 `scanned` (the genuinely image-only PDFs); 1 `broken_text_layer`.

| paper | needs_ocr | gibberish | text_layer_scripts | visual_script | packs selected |
|---|---|---|---|---|---|
| Olfers 1824 | False | 0.302 | {Latin: 1.0} | None | **['cat']** |
| Eschscholtz 1825 | False | 0.191 | {Latin: 1.0} | None | ['deu','deu_latf'] |
| Keferstein & Ehlers 1860 | False | 0.229 | {Latin: 1.0} | None | ['deu','deu_latf'] |
| Chun 1882c | False | 0.258 | {Latin: 1.0} | None | ['deu','deu_latf'] |

These PDFs carry a garbage text layer from some earlier OCR — the dataset README documents this
as a property of the material. corpus trusts it. What lands in `chunks.json` for Olfers 1824:

```
!Oir llt/ne €lHbfl1rr. … S!3Dn ber @r6h eintt1 't.atlbencvd. … ~n bm tropird}tn ~tmn.
```
(= *Die kleine Seeblase. … Von der Größe eines Taubeneies. … In den tropischen Meeren.*)

Three compounding causes:

1. **The gibberish heuristic cannot see blackletter corruption.** `ocr.gibberish_threshold` is
   0.65; these score 0.19–0.30, because the corruption is letter-like Latin characters rather
   than symbol soup. The threshold is effectively unreachable for this failure mode, so
   `needs_ocr` is never set and `deu_latf` is never invoked — while README §Language support
   promises exactly this capability and `scan_detection.json` cheerfully records
   `tesseract_pack_available: true`.
2. **`visual_script` is `null` on all 32 documents.** It sits right next to
   `text_layer_scripts` and is precisely the cross-check that would catch this — page image
   looks like Fraktur, text layer claims clean Latin. Nothing populates it.
3. **Language detection runs on the corrupted text.** Olfers 1824 → `detected_language: "ca"`
   → `tesseract_packs: ['cat']`. So even forcing OCR there would run the Catalan pack.

Measured downstream damage: 3 taxa extracted from Eschscholtz (which gives ~66 chunks to a
systematic treatment of siphonophores), 1 from Olfers; **zero** anatomy lexicon terms in either
though the text names *Kieme*, *Gefäße*, *Herz*; `section_class` null on every chunk. These
papers are invisible to taxon-, lexicon- and topic-driven queries while `corpus status` shows
100% across every stage.

**Proposed fix**
- Populate `visual_script` and gate on disagreement with `text_layer_scripts`; that is the
  robust signal and the plumbing is already half-built.
- Add a Fraktur/blackletter-specific detector for the text-layer path — the `ſ`→`f`, `ch`→`cf)`
  substitution signature is highly regular and cheap to score. Or lower the threshold *for the
  text-layer-trust decision specifically* (distinct from the post-OCR quality gate).
- Run `langdetect` on OCR output rather than a suspect text layer, or refuse a detection whose
  confidence is high but whose language is implausible for the pack set.
- Ship a documented escape hatch: `ocr.force_ocr: true` in config and `corpus run --force-ocr`,
  plus a per-paper override. There is currently none in README, INSTALL.md or the template.
- Add a quality gate that fires when a pre-1950 document takes the `clean_text_layer` path, so
  this at least appears in `corpus status`.

### 1.1a What shipped for 1.1
**Detection.** New `_scanned_page_fraction` asks whether pages carry a single image covering
≥50% of the page — the "is this a scan?" question, independent of the text layer. Above
`ocr.scan_page_fraction_min` (0.40) the document is `scanned` / `raster_page_images` and gets
re-OCR'd. On the reference corpus the two populations are perfectly separated: every scan
scored 0.50–1.00, every born-digital paper exactly 0.00. **29 of 32 papers now OCR, up from 4.**

Page-count fraction, not mean image area: mean area puts a two-page scan with one plate page at
0.50, too close to a digital paper carrying large figures. It samples *across* the document,
not the first N pages — Kawamura 1911a is a born-digital English typescript for pages 0–7 and a
full-page Japanese scan from page 12, and front-only sampling reported 0% raster and skipped
OCR on 13 pages of Japanese.

**Mode.** `--force-ocr` for uniformly-scanned documents; `--redo-ocr` for mixed ones
(scanned fraction 0.40–0.95), which replaces OCR text while leaving genuine digital text alone
so a bound-in typescript isn't rasterized to fix a scanned half. Falls back to `--force-ocr` if
ocrmypdf refuses. Never `--skip-text` — that would preserve the layer we just rejected.

**Language.** A language read off a rejected text layer inherits its unreliability (Olfers 1824,
German Fraktur, was detected as Catalan and routed to the `cat` pack). So the language now comes
from OCRing a 5-page sample:
- Sampled from the middle 15–85% — the front is covers/plates/boilerplate, the tail is
  references, which are a multilingual pile of proper nouns and the worst possible detection
  input. A naive 75% sample lands in the bibliography of a 314-page monograph.
- **Per page, then unioned.** Bilingual volumes are routine here and one verdict per document
  is simply wrong for them. Verified: Kawamura 1911a → `jpn+eng`, Carré 1969 → `eng+fra`,
  Margulis 1976a → `eng+rus`, Gasca & Suárez 1993 → `eng+spa`.
- Tesseract OSD picks each page's script first, then that page is OCR'd with only that script's
  packs. Probing every page with the full 13-pack union cost ~85 s/document for the same answer.
- Pages under 200 characters are skipped; a plate or blank verso shouldn't outvote a body page.

**Two empirical results worth not re-litigating:**
- **Probe DPI must stay 300.** At 200, Kawamura lost its Japanese entirely, Linnaeus 1735 went
  from correctly finding nothing to a confident wrong `ca`, and Bernstein 1934 gained a spurious
  `bg`. The probe is the slowest part of detection and therefore the obvious thing to optimize;
  don't optimize it this way.
- **langdetect ships no Latin profile** — `la` is not among its 55 languages, so Latin can only
  be mis-identified, as a low-confidence Romance language. `probe_language_min_confidence`
  (0.85) routes those to the script-narrowed union, which does include `lat`. That is a safety
  net, not a fix; identifying Latin properly needs a different detector. Relevant because Latin
  diagnoses and pre-Linnaean works are core taxonomic material.

**Discarding wrong-pack output.** Tesseract OSD is not infallible: on Bernstein 1934 it read a
Cyrillic table as Latin, Tesseract transcribed it with Latin packs
(`Ta6auna 4 ... Bron. | NeNe cranguk`), and langdetect called that Catalan at p=0.86 — so
confidence cannot catch this. Gibberish score can: 0.61 for that page against 0.12–0.39 for the
document's genuine Russian and German pages. Pages above `ocr.probe_max_gibberish` (0.50) are
discarded. **The gate applies only to Latin-dominant page text** — `_gibberish_score` counts
≤2-character tokens as suspicious, which is right for Latin-1 mojibake and wrong for CJK, and
gating CJK on it threw away Kawamura's Japanese a second time.

**A correction to the dataset's own documentation:** Bernstein 1934 is authentically bilingual
Russian *and* German (page 18: `Auf Grund der Verteilung dieser Form in der Barents-See kam
Linko zu der Schlussfolgerung...`), a substantial German section of the kind normal for pre-war
Soviet zoology. `siphonophores_smoke/readme.md` describes it as Russian only.

**Smoke set extended** to 35 papers at CWD's suggestion, with bib records and readme updated:
`Kawamura1911a` (digital English front / Japanese scan back — the mixed-volume case),
`Carre1969_Nanomia_tr` (French front / English back), `Hosiaetal2024` (*Nanomia bijuga* figures).

### 1.1b Rebuild results (2026-08-02)
Full rebuild of the 35-paper smoke corpus with OCR enabled: **exit 0, 35 papers / 3,323 chunks
/ 420 figures**, every stage 100%, 31 documents OCR'd (was 4).

Olfers 1824, the same passage, before and after:
```
before:  !Oir llt/ne €lHbfl1rr. … S!3Dn ber @r6h eintt1 't.atlbencvd.
after:   …sehr deutlich und zum Theil abweichend von den früheren Untersuchungen.
         Ueber den sogenannten Giftsporn des männlichen…
```
Keferstein & Ehlers 1860 now OCRs its own title correctly (`Auszug aus den Beobachtungen über
die Siphonophoren von Neapel und Messina angestellt im Winter 185…`, matching the bib record).
Kawamura 1911a's Japanese half reads under `jpn+eng`. Eschscholtz 1825 went from 3 extracted
taxa to 6.

**`zero_references_unexpected` rose 8 → 12, and that is correct, not a regression.** It looked
like force-OCR had destroyed Grobid's reference parse, so I ran Grobid against the original and
the OCR'd PDF for each newly-flagged paper. The references only the *original* layer produced
were fabricated out of mojibake:

| paper | refs from original | refs from OCR'd | what the original's "references" were |
|---|---|---|---|
| De Haan 1827 | 23 | 0 | surnames `Den Nam Wui I .`, `Jdt`, `Ta`, `On` |
| Keferstein & Ehlers 1860 | 12 | 0 | title `mac!}:: ~aber ber. (f5'd)roimmfacf gebilbet ifl…` |
| Chun 1882c | 3 | 0 | — |
| Vanhöffen 1906 | 0 | 0 | no change |
| Beklemishev 1969 | 0 | 0 | no change |

Grobid was hallucinating reference structure from corrupt text and parsing taxonomic names
(`Apolemia Uvaria`, `Abyla pentagona`) as reference titles. Emitting zero is the honest result.
This is also likely relevant to A3 (reconciliation matching nothing): the authority DB was being
fed these.

**One real regression, found and fixed: `_gibberish_score` is Latin-centric.** Yamamori 2014
OCR'd correctly under `jpn` and was still flagged `gibberish_after_ocr` at 0.55, because the
"token of ≤2 characters is suspicious" rule is meaningless for CJK. A document-level
Latin-dominance test did *not* fix it — that paper is 88% Latin by character (English abstract,
references, taxonomic names) with Japanese body text. Fixed at the root: `_gibberish_score` now
excludes CJK tokens from both numerator and denominator. Yamamori drops to 0.45; the plate-only
Quoy & Gaimard 1834 still correctly flags at 0.58.

**This assumption has now bitten three times in one session** — the OCR language probe's
gibberish gate, the `gibberish_after_ocr` quality gate, and the original scan-detection
threshold that could not see blackletter corruption. Worth treating as systemic: any heuristic
tuned on Latin prose should be audited before it is applied to a corpus that is deliberately
multi-script. A residual instance: Chen et al. 2015 (Chinese) scores 0.75 even with CJK tokens
excluded, because the non-CJK remainder is mostly short fragments and numerals. It is harmless
today only because that paper is born-digital and the gate is gated on `needs_ocr`.

**Still open from 1.1:** the `visual_script` field remains unpopulated on the non-raster paths,
and there is still no quality gate that fires when a pre-1950 document is trusted as
`clean_text_layer`. Both are now much less load-bearing, since the raster check catches the
population they were meant to catch.

### 1.2 `corpus prefetch` missed a model, and the miss degraded chunking silently — **FIXED**
`corpus prefetch` managed three models and asserted *"every required model is already cached;
nothing to do"*, while the pipeline also loads
`sentence-transformers/all-MiniLM-L6-v2` (docling's `HybridChunker` tokenizer) at chunking time.

Verified by removing it from the cache and following the documented offline recipe
(INSTALL.md §Offline hosts — prefetch where there is network, `HF_HUB_OFFLINE=1` where there
isn't):

```
WARNING pipeline.chunking: HybridChunker failed (We couldn't connect to
  'https://huggingface.co' ...); falling back to naive char chunker
INFO pipeline.chunking: Naive chunker produced 1 chunks
✓ extract completed in 10.7s
═══ All 1 step(s) succeeded in 10.7s ═══
```

Same paper, same config: **16 chunks with the tokenizer, 1 chunk without.** Exit 0 either way.
On a cluster this degrades retrieval across the entire corpus with nothing but a per-paper
WARNING to show for it — the same failure shape as #139.

Fixed on this branch:
- `pipeline/prefetch.py` — `prefetch_docling` now chunks the sample document as well as
  converting it, so whatever tokenizer the pinned docling wants gets cached; the tokenizer repo
  is added to `DOCLING_REPOS` so `corpus check` counts it and `corpus prefetch` stops
  short-circuiting on an incomplete set.
- `pipeline/chunking.py` — the fallback now logs at ERROR and names the consequence and the
  remedy.

Verified: `corpus check` now reports `1 of 4 not in the local HuggingFace cache
(sentence-transformers/all-MiniLM-L6-v2)`, and `corpus prefetch` fetches it.

**Still open:** the naive-chunker fallback should surface in `corpus status` as a quality flag,
not only in per-paper logs. `chunks.json` already records `chunker_name`, so the data is there.

### 1.3 `install_tessdata.sh` left low-accuracy OCR models in place — **FIXED**
README, INSTALL.md and environment.yaml all asserted "the conda-forge `tesseract` package ships
only the English LSTM model". Not true any more: conda-forge `tesseract` 5.5.2 bundles 158
`.traineddata` files, and they are **tessdata_fast** builds — verified byte-exact against
upstream `tesseract-ocr/tessdata_fast` (`deu` 1,525,436 B; `rus` 3,861,738 B) versus
tessdata_best (`deu` 8,628,461 B; `rus` 15,301,764 B).

Because the installer skipped any file already on disk, **a fresh conda env downloaded exactly
one pack** — `deu_latf`, the only code conda-forge doesn't ship — and printed "already present,
skipping" for the other twelve, under a banner reading "Source: tessdata_best (high-accuracy
LSTM)". `eng` was not in the script's language list at all, so the pack README calls the
always-appended final fallback was never upgraded either. `corpus check` cannot see any of
this; it reports "all 13 configured packs installed".

Fixed: the installer records what it fetched in
`$CONDA_PREFIX/share/tessdata/.corpus_tessdata_best` and replaces any pack not in that marker
("replacing non-best copy"), preserving idempotence without trusting mere file presence; added
`--force`, `--help` and unknown-option rejection; added `eng`; corrected the three docs.

Note this is currently **latent** — because of 1.1, almost nothing in this corpus actually
reaches Tesseract. Fixing 1.1 is what makes 1.3 matter.

---

## P2 — first-run experience

### 2.1 `corpus run --dry-run` on a fresh corpuscle printed 6 ERRORs — **FIXED**
Seconds after `corpus check` said *"ready: `corpus run` should succeed on this host"*,
`--dry-run` emitted six ERROR lines and six warnings. Five independent guards treated
"corpuscle not built yet" as failure. Made each dry-run-aware so they report the real plan
("would embed 0 hash(es)", "0 hash dirs"): `pipeline/main.py` (the #139 taxonomy precondition —
a dry-run writes nothing, so it cannot produce the empty `taxa.json` the guard exists to
prevent; real runs unchanged), `pipeline/embed.py`, `pipeline/intext_citations.py`,
`pipeline/taxon_mentions.py`, `bib/authority.py`, `bib/reconcile.py`, plus an upfront banner in
`pipeline/orchestrator.py`.

Before: 6 ERROR + 6 WARNING. After: `═══ All 7 step(s) succeeded in 2.2s ═══`.

### 2.2 `--dry-run` wrote to disk — **FIXED**
It created `output/`, `output/documents/` and `output/vector_db/` while printing "No files
written". Beyond the broken promise, this made the *second* dry-run behave differently from the
first, since the guards in 2.1 test for those directories. `pipeline/main.py` now computes the
paths without creating them under `--dry-run`. Verified: only `config.yaml` remains afterwards.

### 2.3 `~` in config paths was not expanded — **FIXED**
`input_pdfs: ~/data/...` produced
`[fail] input_pdfs: /home/claude/corpora/siph_smoke/~/data/... not found` — the tilde glued on
as a relative path component. `_resolve_against` in `pipeline/cli.py` now expands it (single
choke point, covers every path field); the template documents that `~` works and `$VARS` don't.

### 2.4 The config template's taxonomy block was self-contradictory — **FIXED**
It said "uncomment EVERY line below", over a block that would then yield `source: worms` *plus*
a `path:` — an invalid mix, defaulting to the network-walking source the README elsewhere tells
you to avoid. Replaced with two clearly-alternative blocks, `dwca` first as the recommended
default.

### 2.5 `corpus status --filter-gate <name>` silently ignores the filter — **FIXED (hint only)**
The report prints `List affected papers with: corpus status --filter-gate <name>`; running
exactly that reprints the whole report, because the flag only takes effect alongside
`--list-hashes`. Fixed the hint. **Still open:** better to make `--filter-gate` on its own
list the affected papers, since that is what the text promises.

### 2.6 Linux `docker` group not mentioned — **FIXED (doc)**
`apt install docker.io` leaves the user unable to reach the daemon; every `docker compose` line
in the README then fails with `permission denied ... /var/run/docker.sock` — after they have
already installed an 8 GB env. Added the `usermod -aG docker` step and a `hello-world` check.

### 2.7 Disk budget omitted the env; prefetch implied offline-readiness — **FIXED (doc)**
README §Computational requirements counted ~5 GB of models but not the **8.1 GB** conda env.
Now "~15 GB before a single PDF is processed". Also noted that prefetch does not make a run
network-free — HF revision checks still go out; `HF_HUB_OFFLINE=1` is the setting that pins it.

---

## P3 — MCP tool surface (found by exercising the built corpuscle)

### 3.1 The README's own example figure query returns nothing relevant
README §Example uses advertises *"Give me every figure showing nectophore morphology in*
Nanomia bijuga*."* In practice `get_figures_for_lexicon_term(anatomy, "nectophore")` returns
~40 figures and **not one is *Nanomia***, because it matches captions only. The best actual hit
— Totton & Bargmann text-fig. 37, whose body text reads "distinguished from *N. bijuga* by its
nectophores (text-fig. 37)" — has a caption of just "*Nanomia cara* (A. Agassiz)". A client
that trusts the tool concludes the corpus has no *Nanomia* nectophore figures.
**Fix:** let figure search consider `figure_refs` from body chunks, not captions alone; or
document the caption-only semantics prominently in `MCP_TOOLS.md` and soften the README example.

### 3.2 Caption panel-parsing invents panels from author initials and credit lines
`panels_from_caption` with `kind: "period"` splits on `<letter>.` without excluding name
abbreviations. **17 such records** in this corpus:
- `{'label': 'C', 'description': 'Munro'}` ← "Photo credit to C. Munro"
- `{'label': 'A', 'description': 'Agassiz'}` ← the authority "(A. Agassiz)"

These become bogus ROIs and mislead `list_figure_rois` consumers.
**Fix:** exclude single-initial-plus-surname patterns and trailing credit/authority clauses.

### 3.3 No language field anywhere in the tool surface
`scan_detection.json` records `detected_language` per paper; nothing surfaces it — not
`list_papers`, `get_papers`, nor `corpus_summary`. Asked what languages the corpus contains, a
client must infer from titles. For a project built around multilingual historical literature
that advertises cross-lingual retrieval, this is a conspicuous gap, and it blocks the obvious
"show me only the Russian papers".
**Fix:** add `language` to the paper metadata the bundle carries and expose it in `list_papers`
/ `get_papers` / `corpus_summary`, with a `language=` filter.

### 3.4 Figure licensing fields are null
`get_figure` returned `license: null` / `license_url: null` for every figure inspected, with
only `attribution` populated, despite `licensing.pd_cutoff_years: 95` and a corpus that is
mostly pre-1930. Worth confirming whether the strict `manuscript` / `presentation` profiles can
ever clear a figure in practice.

---

## P1.4 — OCR silently dropped whole pages to a per-page timeout — **FIXED**
Found while verifying the OCR rebuild, and the most serious defect in the review. OCR was
losing entire pages, **nondeterministically**, and reporting success. Linnaeus 1735, identical
command and input across two runs:

```
run A  [43, 414, 5027, 6420, 10575, 3796, 6012, 0, 9758, 8462, 6124, 0, 0]
run B  [43, 414, 5027,    0,     0, 3796, 6012, 0,    0,    0,    0, 0, 0]
```
56,631 characters became 15,292. Both logged `OCR completed successfully`, both exited 0, and
`corpus status` reported 100% on every stage.

Reproducing with stderr captured named it immediately:
```
12 [tesseract] took too long to OCR - skipping
13 [tesseract] took too long to OCR - skipping
 9 [tesseract] lots of diacritics - possibly poor OCR
   Suppressing OCR output text with improbable aspect ratio
```
`--tesseract-timeout`, documented as "give up on OCR after the timeout, but copy the
preprocessed page into the final output" — a blank page and a clean exit. The pipeline never
set it, so it inherited ocrmypdf's default, far too tight for a dense 300-dpi scan with seven
language packs loaded. Load-dependent, so it does not reproduce reliably.

Fixed three ways: the timeout is set explicitly (`ocr.tesseract_page_timeout`, 900 s/page);
ocrmypdf's stderr is logged **on success**, which is the single line that made this invisible;
and a post-OCR check names pages that came back with no text.

With the timeout raised, Linnaeus yields **92,616 characters across all 13 pages** — meaning
the *better* of the two runs above was still missing 39% of the document, and pages 12–13 were
blank in both. This was never a two-run fluke; it was constant, unmeasured loss.

Exposure scales with the re-OCR change: 31 documents take this path where 4 did before. On a
cluster with array jobs competing for cores it would be worse and equally silent.

## Final verification (2026-08-02, 35-paper rebuild with every fix)

Clean rebuild, 111 min, exit 0, all stages 100%, no failures.

| check | before | after |
|---|---|---|
| corpus text extracted | 1,943,374 chars | **2,206,893** (+13.6%) |
| Linnaeus 1735 alone | 38,349 | **301,833** |
| bogus author-initial panels | 3 | **0** |
| chunks with a `section_class` | 357/3211 (11.1%) | **555/3412 (16.3%)** |
| `section_class: description` | 18 | **193** |
| `section_class: methods` | 1 | **23** |
| `gibberish_after_ocr` flags | 2 (one false) | **1** (the genuine plate-only volume) |
| papers OCR'd | 4 (pre-review) | **31/35** |

`corpus serve --check` green. `scan_detection.json` reaches all 35 distilled bundles, so the
language surface works end to end: `by_language` = en 16, de 7, fr 7, ru 5, sv 2, it 2, ja 2,
es 1, nl 1, pt 1, and `language="ru"` returns exactly the five Russian papers the dataset
documents.

Two defects found *during* this verification and fixed:

* **The CJK gibberish fix was insufficient.** Yamamori 2014 still flagged, because the quality
  gate scores `text.json` rather than chunk text, and excluding CJK *tokens* is not enough when
  the non-CJK remainder is `['##', 'Li', '\_', "『'", '=']`. On a document 55–63% CJK the
  heuristic has no validity in either direction, so `_gibberish_score` now returns 0.0 above
  `_CJK_SHARE_UNSCORABLE` (0.30) instead of a misleading number. Safe only because
  `_scanned_page_fraction` now answers "is this a scan?" independently — this heuristic is no
  longer the sole guard.
* **A bug in the new language surface.** `_scan_facts` fell back to `detected_language` even
  when `language_trusted` was False, laundering a value the pipeline had explicitly rejected
  into the served API — Linnaeus 1735 was reported as Catalan. It now reports no language.

**Known caveat, not fixed:** langdetect has no Latin profile, so Latin text surfaces as a
low-confidence Romance language. Tilesius 1814 (Swedish, with Latin passages) reports
`ca` among its languages. Consumers of the `language` field should know that `ca` on a
pre-modern work often means Latin.

**The Latin-centric assumption in `_gibberish_score` bit five separate times in this review** —
scan detection's blackletter blindness, the OCR language probe's page gate, the
`gibberish_after_ocr` quality gate, the CJK token exclusion, and finally the text.json scoring
path. Any heuristic tuned on Latin prose needs an explicit decision about what it does on a
corpus that is deliberately multi-script.

## T1/T2 integration tier — run for the first time, with a control

The `corpus_required` tier had never been run on this branch (T0 deselects it). Run against a
freshly rebuilt demo it gave 5 failures. A control — demo rebuilt at the branch point
`7a20d80`, my changes absent — gave 5 failures too, but not the same 5:

| test | branch point | this branch |
|---|---|---|
| `test_title_appears_in_text[af043530e5dd]` (Marrus) | FAIL | FAIL |
| `test_title_appears_in_text[dde93d15a5e8]` (Stepanjants) | FAIL | FAIL |
| `test_first_author_surname_in_text[dde93d15a5e8]` | FAIL | FAIL |
| `test_references_match_corpus_papers[dde93d15a5e8]` | FAIL | FAIL |
| `test_references_match_corpus_papers[4fe914163f59]` (Pugh) | **FAIL** | **passes** |
| `test_title_appears_in_text[ef8482d9cb44]` (Schneider) | passes | **FAIL** |

Four pre-existing; **one fixed, one introduced**. Both of the movers are informative.

**The introduced failure is the OCR change working correctly on a mis-classified document.**
`demo/Schneider1891.pdf` has a full-page image on *every* page plus embedded fonts — a scan
carrying a text layer, not born-digital. The README describes it as "born-digital with a clean
text layer" and builds its explanation of what the default demo run exercises on that claim.
The claim is wrong, and the "clean" layer is mediocre third-party OCR:

```
control (trusted text layer)      branch (corpus re-OCR)
Prof. J. Victor (Jarus       ->   Prof. J. Victor Carus        ✓
Verlag von "rilhclm Engelmann ->  Verlag von Wilhelm Engelmann ✓
No. 353-:180.                ->   No. 353- 380.                ✓
1091.                        ->   1591.            (both wrong)
5,696 chars                  ->   6,009 chars
```
corpus's own OCR recovers the author, publisher and issue numbers. It also introduces a few
noise lines (`eee`, `AE`, `IO`) from decorative rules, and still misreads the year. Net
improvement.

`test_title_appears_in_text` fails anyway, because it exact-substring-matches a BibTeX title
against docling's extracted text. **Any** re-OCR breaks it regardless of whether quality rose
or fell — it tests string identity, not correctness. The same brittleness explains two of the
four pre-existing failures: Marrus fails on a single space before a comma
(`marrus claudanielis, a` vs `marrus claudanielis , a`), with the title plainly present.

**Recommended, not done here** (test-design changes deserve your call): normalize whitespace and
punctuation before comparison, or score similarity rather than requiring a substring. As
written the test will fail intermittently on any docling upgrade.

**Also fix the README's demo description** — Schneider 1891 is a scan, so the sentence claiming
a default demo run does not exercise OCR on it is no longer true.

## P2.8 — lexicon translations match only uninflected forms  [PLAN — verified]
Re-exercising the rebuilt corpus over MCP surfaced that Olfers 1824 and Eschscholtz 1825 still
extract **zero** anatomy terms, while Vanhöffen 1906 — also German — extracts eight. So German
matching works in general; something narrower is wrong.

The MCP client's diagnosis ("the lexicon has no German surface forms, the #143 synonym
resolution covers English↔Latin and nothing else") is **wrong**: `demo/lexicon.yaml` carries 17
German forms plus French and Russian, and the pipeline consumes them. The real cause is one
step further in:

| text | matches? | why |
|---|---|---|
| `Luftblase` | yes | listed under `translations.de` |
| **`Luftblasen`** — as actually printed in Eschscholtz | **no** | inflected form not listed |
| `pneumatophore` | yes | canonical |
| `pneumatophores` | yes | English plural hand-listed in `synonyms` |

`_build_lexicon_matcher` matches whole words against an enumerated surface-form set. English
plurals only work because someone typed them into `synonyms`. The `de` / `fr` / `ru`
translations list base forms only, so in German — which inflects heavily — most real
occurrences are missed. Eschscholtz's actual text reads `Luftblasen`, and the paper's other
anatomical vocabulary (`Saugmägen`, `Fangfäden`, `Schwimm- und Athmungshöhle`) is not in the
lexicon in any form.

The effect is worse than a low score: a German paper reports anatomy coverage of exactly zero,
which reads as "nothing to see here" rather than "not indexed".

**Not fixed** — two defensible routes and the choice is a design decision:
1. Enumerate inflected forms in the lexicon (explicit, no false positives, tedious, and each
   new language repeats the work).
2. Give the matcher light per-language suffix tolerance for non-English surface forms — real
   false-positive risk, so it wants test coverage before it ships.

Whichever, the lexicon's non-English vocabulary is also just thin for pre-1900 German, which is
a content problem independent of matching.

**Also from that session, worth checking separately:** species-level taxa are not extracted from
abbreviated binomials (`Ph. producta`), so Olfers 1824 — a paper whose entire purpose is
erecting *Physalia producta* — yields no species-level taxon. I did not verify this one.

**And one claim from that session that is _not_ a defect:** `search_taxon("Porpita")` and
`search_taxon("Velella")` return `not_found`. Both are Porpitidae (Anthoathecata), correctly
outside a Siphonophorae-rooted Darwin Core snapshot — `demo/instructions.md` exists partly to
tell the model exactly this. Absence is correct scope, not missing coverage.

## P3 corrections — two findings withdrawn after measurement

**A3 (reconciliation) is not broken.** The claim that bibliographic resolution was "effectively
non-functional" came from an MCP client's inference over one paper's reference list, and I
repeated it without measuring. The authority DB over 710 citations: `doi_exact` 233,
`title_fuzzy` 7, `alias_exact` 12, `new_work` 458. Fuzzy title matching exists and fires.
The `matched 0 / no_candidates 27` line is a *different* step (merging ghost cited-references
onto corpus papers), and in a 35-paper corpus spanning 1594–2026 most papers genuinely are not
cited by the others.

The 5 `low_score` near-misses are mostly **correct rejections** — same author and year,
different works:

| corpus paper | best ghost candidate | verdict |
|---|---|---|
| Kawamura 1911 | `shidarezakura kurage and nagayoraku kurage` | different paper |
| Carré 1969 | `etude du developpment larvaire de sphaeronectes` | different paper (corpus one is *Rosacea villafrancae* sp. n.) |
| Chun 1882 | `ueber die cyclische entwickelung…` | different paper |
| Totton 1965 | *(empty title, score 0)* | unverifiable |
| Vanhöffen 1906 | `siphonophoren nordisches plankton` | plausibly the same work |

Lowering the threshold would merge distinct works, which is worse than leaving them separate.
**No change made.**

**3.4 (licensing) is not a bug.** With `pd_cutoff_years: 95` the cutoff is 1931. All 16
pre-1931 corpus papers are correctly `public_domain` via `license_source: age_based_pd`; the 17
`no_record` ones are all 1932 or later, where conservative "unknown" is right because they may
still be in copyright. `license: null` means "no explicit license string" — the determination
lives in `publishable` / `license_source` / `clearance_state`, surfaced under a strict profile
or `include_licensing=True` (#154). Worth a line in MCP_TOOLS.md, not a code change.

**A2 (`section_class`) was also mis-framed** — see §1.1b. The 89% null rate is mostly correct
and the multilingual patterns work; only a handful of genuine misses needed fixing.

The pattern across all four: an MCP client's plausible-sounding inference is not evidence.
Every one of these needed measuring against the artifacts before acting, and three of the four
would have been actively harmful to "fix".

## P4 — polish

| # | Finding |
|---|---|
| 4.1 | **README and `corpus --cite` gave different citations** — different first author *and* different title, same DOI. **RESOLVED by CWD:** the citation is now "Church, S. H., Mańko, M. K., Zapata, F., & Dunn, C. W. (2026). Extracting AI agent-accessible data from biodiversity literature with corpus." Applied to `CITATION.cff`, its packaged copy `pipeline/CITATION.cff`, and README §Citing corpus; `corpus --cite`, `corpus --cite=bibtex` and `bundle_manifest.json` now agree byte-for-byte. Manko → Mańko throughout. **Follow-up:** `pipeline/CITATION.cff` is a hand-maintained byte copy of the root file with nothing enforcing that they match — that is how these drifted apart. A test asserting the two are identical is cheap insurance and is not yet written. Bundles built before this change still carry the old citation; it is stamped at build time. |
| 4.2 | Third-party warnings front-run every command: torch's `CUDA initialization ... Error 804` prints twice before any corpus output on a driver-mismatched host, and `[transformers] Token indices sequence length is longer than ...` leaks to the console during chunking without a timestamp or module prefix. Both are benign; both look like errors. Suppress or route through the logger. |
| 4.3 | Long per-document stages give no progress signal — the log went silent for 3m20s during docling layout on a 27-page scan, far longer on the 314-page monograph. A periodic heartbeat would distinguish "working" from "hung". |
| 4.4 | README §Instructions says `<corpuscle>/instructions.md` while §Corpuscle layout defines `<corpuscle>/` as the tree holding `documents/`, i.e. `output_dir`. It is actually read from the config directory. Pin down one term. |
| 4.5 | `dev_docs/clean_install_walkthrough.sh` §5 puts `corpus run --dry-run` immediately after `corpus check`, which is what made 2.1 so visible. Worth re-running the walkthrough end to end after these fixes. |

---

## Suggested sequencing

1. **1.1** — the only finding that makes the shipped corpus substantively wrong. Everything
   else is either fixed or cosmetic by comparison. Ship-blocking in my judgement.
2. ~~**4.1**~~ — done; see the table.
3. **3.1 / 3.3** — the tool surface not matching the README's advertised queries.
4. **1.2 residual** (chunker fallback as a status flag), **2.5 residual**, **3.2**.
5. **P4** polish.

## Changes already on this branch

```
INSTALL.md                    README.md                  environment.yaml
bib/authority.py              bib/reconcile.py           pipeline/chunking.py
pipeline/cli.py               pipeline/config.template.yaml
pipeline/embed.py             pipeline/intext_citations.py
pipeline/main.py              pipeline/orchestrator.py   pipeline/prefetch.py
pipeline/status.py            pipeline/taxon_mentions.py tools/install_tessdata.sh
```

All 696 T0 unit tests pass (`pytest -m "not corpus_required"`). No new pyflakes findings.

---

## Addendum — findings from a second round of MCP prompts

Four more, all reproduced against the built corpuscle.

### A1 — `get_chunks_for_topic` documents `score` as a similarity but returns a distance — **FIXED (doc)**
`mcpsrv/tools/chunks.py:131` described the tool as "LanceDB cosine **similarity**", while
line 218 assigns `r.get("_distance")` to the field named `score` (the code comment there even
says "LanceDB returns cosine distance"). Observed behaviour: rows sorted **ascending**
(0.723 → 0.911) with values **exceeding 1.0** (up to 1.124), impossible for cosine similarity,
normal for cosine distance.

That docstring is the tool description every MCP client reads. A client that sorts descending,
or thresholds on `score > 0.8`, gets exactly inverted results — silently, with plausible-looking
output. Highest blast radius of anything in the tool surface.

Fixed the description to state that `score` is a distance, lower is closer, values run 0–2, and
rows arrive pre-sorted. Renaming the field would be a breaking API change, so it stays `score`.

### A2 — `section_class` is null on 90.6% of chunks  [PLAN]
Measured across all 3,272 chunks:

| section_class | chunks | share |
|---|---|---|
| **None** | **2,964** | **90.6%** |
| references | 153 | 4.7% |
| abstract | 69 | 2.1% |
| introduction | 28 | 0.9% |
| conclusion / description / discussion / acknowledgements | 55 | 1.7% |

Meanwhile **98.2% of chunks carry non-empty `headings`** — so the input is there and
`classify_section` simply isn't matching it. The classified minority is almost entirely English,
which points at an English-only heading vocabulary on a corpus that is deliberately 13
languages. `get_chunks_by_section` is close to unusable as a result, and section context is
missing from every search result.

### A3 — bibliographic reconciliation matches almost nothing  [PLAN]
The build's own reconcile step, from `run.log`:

```
═══ Reconciliation complete ═══
  matched        0
  ambiguous      0
  low_score      4
  no_candidates  25
  missing_text   0
```

**Zero matches across 32 papers.** The MCP session corroborated this from the other side: a
155-reference list from Totton & Bargmann 1965 resolved 100% as `new_work`, so
`get_missing_references` reports works that *are* in the corpus as missing, and Totton &
Bargmann 1965 is itself split into three authority rows
(`corpus:totton|1965|a synopsis of the siphonophora`,
`corpus:bargmann|1965|british museum natural history`, and an empty `corpus:totton|1965|` stub).

The one edge that did resolve did so via DOI (`alias_exact`, Mańko 2020 → Totton 1965), which
suggests exact-key matching works and everything pre-DOI does not. For a corpus that is
overwhelmingly pre-DOI OCR'd literature, that is the wrong default strategy: it needs fuzzy
title matching plus author-surname normalization (Cyrillic homoglyph folding, umlaut
transliteration). `rapidfuzz` is already a declared dependency for exactly this.

Compounding it, Grobid collapsed Totton's ditto-mark reference style (`--1879.`, `--1897b.`)
so single records fuse two or three references (`b25` = Chun 1882 *and* 1885; `b124` = Russell
1877 *and* Sars 1846 *and* 1857). Some entries are unrecoverable from that parse.

This is the root cause of the `get_missing_references` output being unusable as shipped.

### A4 — `resolve_reference` fails on any multi-author query  [PLAN]
`"Totton Bargmann 1965 Synopsis of the Siphonophora"` → `not_found`, with
`parsed_author: "Totton Bargmann"`. Same for `"Keferstein Ehlers 1860 ..."` and
`"Quoy Gaimard 1834 ..."`. Single-author queries resolve fine. Everything before the year is
taken as one surname. Typing both surnames is the natural way to look up a two-author work, so
this fails on first use.

### A5 — `get_citation_graph` default drops two-thirds of a hub work's edges  [PLAN]
Defaults truncate at 50 edges with `truncated: true`, cutting Totton's 155 references off
alphabetically at "Jacobs 1937". The flag is honest and `MCP_TOOLS.md` calls the defaults
"generous", but alphabetical truncation on the corpus's central monograph is a poor default —
rank by `cited_by_count` (which the tool already does for per-node fan-out) rather than by key
order.

### What worked well
Worth recording, since the list above is all problems:
- **Cross-lingual retrieval genuinely works.** English `"pneumatophore structure"` returned a
  Russian Stepanjants 2014 passage at rank 4; a German query put Vanhöffen 1906 and Chun 1898
  at ranks 1–2 while still retrieving English Totton; an English paraphrase avoiding the
  technical term surfaced German and Italian passages. BGE-M3 is doing what the README claims.
- **Taxonomy and synonymy resolution are correct.** *Apolemia* resolved to five accepted
  species with accurate per-species paper coverage and mention counts, and the three uncovered
  species were correctly explained by their post-1965 description dates.
- **`get_figure_image` works** — returned a real PNG that rendered inline.
- **`format_citations` behaves**, including its provenance warnings on reconciled entries.
- **Retrieval crowding is a real effect worth documenting:** 19 of 20 results for
  "pneumatophore structure" came from Totton & Bargmann 1965, which holds 31% of corpus chunks;
  Chun 1898 — a paper titled *Über den Excretionsporus an der Pneumatophore von Physophora* —
  missed the top 20 by 0.02. `get_lexicon_term_dossier` correctly lists all 7 papers and is the
  better entry point for coverage questions. Worth a line in `MCP_TOOLS.md`.

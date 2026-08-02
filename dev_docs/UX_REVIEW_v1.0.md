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

### 1.1 OCR is skipped on 29 of 32 papers, including all four Fraktur ones
The headline multilingual-OCR capability never runs on the material it exists for.

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

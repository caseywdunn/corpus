# Source group boundaries (#317)

Two source cases were found by applying the issue's exact repeated-author
detector to all 35 completed gold TEIs and the four complete source preparations
retained for #314. They are regression cases selected after seeing failures,
not an independent holdout or a new fixture corpus.

- **Totton1965a**, physical page 51: Grobid puts a yearless author-list ref at
  y191 beside refs at y419, with only whitespace between them in TEI. Recovering
  their combined source interval correctly refuses the intervening unrelated
  prose. Separating distant source regions recovers the printed
  `Leuckart (1853, 1854), Gegenbaur (1853), Vogt (1854)` without deleting authors
  globally or attempting to reconstruct the prose Grobid omitted.
- **Tung2003**, physical pages 17 and 19: Grobid joins refs across an intervening
  page. Its page17 boxes also overlap adjacent printed glyphs, producing the
  source interval `。Margulis（1980, 1984）綜`. Each extra boundary glyph straddles
  the corresponding outer box edge. The unique complete parenthesized source
  marker has exactly the raw observations' ordered author/year signature.
  The recovered marker retains fullwidth source parentheses. Its receipt keeps
  the original interval and both excluded boundary glyphs. Fully contained
  neighboring text, initial/suffix letters, digits, author/year disagreement
  and multiple markers are conservative refusal controls.

Each XML is an actual contiguous run of reference elements and their intervening
text tails from the saved full TEI; only the final tail is omitted. Matching
bibliography subtrees are unmodified. Each JSON pins original and prepared PDF
hashes, full/fragment TEI hashes, original metadata, preparation/Grobid receipt,
and exact PDF characters/bounding boxes in the containing intervals. PNGs are
small prepared-source crops reviewed against those character captures. No full
PDF/page or invented source wording is included. The new tests use the actual
coordinate reader against captured glyphs, actual reference parsing and phase1/2
authority, then the production in-text/excerpt functions. Existing Fraser,
Mapstone and Oderberg regressions still cover the named original source cases.

Text evidence and work identity remain distinct. Tung's actual reference parser
splits `R. Ya.` into a second author, so the Margulis 1980/1984 links remain
explicitly unresolved. Totton's reference list also contains unrelated upstream
parse defects. Tests and replay preserve all raw observations and do not invent
links or claim that repairing these markers adjudicates their bibliographies.

## Complete saved-artifact replay

The compact receipt `replay_2026_09_22.json` records **39 documents, 886 citation
paragraphs, 918 immutable reference observations**. The exact issue detector
falls from **22 baseline paragraphs to zero**. Before this follow-up, 20 of those
22 were already repaired; the two remaining hits above were source-reviewed and
fixed. All 22 in-text paragraphs are replayed. **16** also appear through actual
resolved authority links in `get_excerpts_citing`; six have no eligible resolved
link and are retained as such. No artificial edge was inserted. Authority
materialization maps 913 observations; five are withheld by its existing quality
policy. Repeating phase1/2 changes no SQLite content.

Reproduce with the recorded saved input roots (or byte-identical copies):

```sh
python -m tools.qc.citation_paragraphs \
  --input /tmp/corpus-v15-acceptance/gold/output \
  --input /tmp/corpus-issue314-scanned-preparation/build \
  --input /tmp/corpus-issue314-full-preparation/build \
  --output /tmp/citation-replay-new
```

The output path must be new and separate from the inputs. The harness verifies
TEI/PDF receipts and equality of freshly parsed and saved complete reference
lists; it creates only a new authority/artifact tree. Every original observation,
regenerated paragraph and full measurement remains in that output. The committed
receipt keeps source hashes and per-paragraph before/after hashes instead of
copying all measured prose. The historical scratch paths are provenance, not
portable dependencies of the committed regression tests.

The first reporting attempt used a nonexistent `has_more` response key after
successful materialization; the harness now follows the actual `next_offset`
contract. The successful receipt comes from a separate second output directory;
no original OCR, Grobid, or saved build was rerun or modified.

These are saved-artifact measurements using current code, with actual prior
upstream producers named in the receipt. They are not fresh full-corpus
extraction, a replacement September deployment, live MCP transport, or a score
for general text correctness. The detector misses other paragraph defects. The
metadata policy fingerprint invalidates paragraph materialization while allowing
unchanged OCR/extraction and a matching TEI cache to be reused.

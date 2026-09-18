# Citation source evidence (#309, #317)

These are small paragraph/reference fragments and PDF character coordinates,
not another corpuscle. They come from the three source PDFs named by the audit:

| Source | PDF SHA-256 prefix | Physical page(s) | Regression |
| --- | --- | --- | --- |
| Pugh1974.pdf | `b5a7af6140ca` | 5 | Fraser (1961, 1967), with separate reference targets |
| Oderberg2020.pdf | `c59691fa845a` | 9–10 | Totton and Mackie 1960 / Bardi and Marques 2007, incorrectly split and linked to the 1965 Synopsis |
| Mapstone2009.pdf | `4efbb4af134b` | 12 | Shared-author Pugh years, letter lists, and the 1992a–c range |

Each XML file contains the original paragraph subtree and relevant reference
subtrees from a fresh **lfoppiano/grobid:0.8.1** run on the source PDF, with
`consolidateHeader=0`, `consolidateCitations=0`, `includeRawCitations=1`, and
`teiCoordinates=p,ref`. The JSON records the full source PDF hash, full captured
TEI hash, fragment hash, and source PDF characters/bounding boxes, limited to
citation-bearing intervals on the indicated pages. Coordinates are rounded to
three decimal places. The regression invokes the actual coordinate-to-text
implementation against those captured characters; a separate generated-PDF
integration test exercises the PDF reader itself.

The fresh TEI reproduces the two audit mechanisms directly. It inserts
`Fraser ( , 1967) )` into the second Fraser reference, despite the source's clean
group, and splits the Oderberg group into `(Totton and`, `Mackie 1960, Bardi and`,
and `Marques 2007)`. The first clipped piece targets `b30` (1965). The source
contains the complete Totton and Mackie 1960 span, whose reference is `b31`.
These are observed upstream outputs, not invented TEI test strings. Grobid's
[upstream issue #830](https://github.com/grobidOrg/grobid/issues/830) also documents
author expansion in TEI reference text.

Source repair reconstructs each citation group from the prepared PDF while
keeping literal TEI surfaces, targets, and coordinates in `tei_observations`.
Separate author/year spans allow several targets to share one printed author
or letter range. Complete surname/year matching is conservative: unresolved
and ambiguous references retain their text and candidate evidence instead of
receiving a guessed target. In the captured Mapstone reference list, both
`b20` and `b521` describe Pugh1999b; the parser reports that ambiguity. Other
upstream reference defects also remain visible: the Pugh1974 bibliography gives
Leloup1955 a structured year of 1910, and some Totton entries are merged. This
repair does not claim to resolve those bibliography identities.

## Source replay observations

The audit's repeated-author signature was replayed over **all citation-bearing
paragraphs of these three fresh TEI captures**, before and after source repair:

| Source | Citation paragraphs | Signature hits before | After |
| --- | ---: | ---: | ---: |
| Pugh1974 | 101 | 25 | 0 |
| Oderberg2020 | 25 | 0 | 0 |
| Mapstone2009 | 645 | 20 | 0 |

There were no remaining signature hits in these captures to inspect manually.
These are signature counts on these specific observations, not error rates or
proof of corpus-wide repair. Paragraphs without matching PDF evidence remain
explicitly TEI-derived. The target-validation statuses after replay were:

| Source | Validated author/year | Unresolved author/year | Ambiguous author/year | Unverified TEI |
| --- | ---: | ---: | ---: | ---: |
| Pugh1974 | 175 | 138 | 2 | 19 |
| Oderberg2020 | 41 | 8 | 0 | 3 |
| Mapstone2009 | 2,102 | 1,198 | 384 | 3,089 |

Unresolved/ambiguous counts are **review candidates**, not verified target
errors. Raw targets can be withheld because bibliography metadata is wrong,
missing, or duplicated. In Oderberg, the only previously resolved marker whose
old target is removed is the demonstrated clipped `b30`; complete 1960 and
1965 occurrences remain linked. The integration regression checks the rebuilt
artifact through citation excerpts, in-text pagination, and citation graph
queries, including independent valid citations of both works.

The normal metadata stage now requests ref coordinates and binds the TEI cache
to that request policy. A release rebuild must regenerate legacy TEI to obtain
those coordinates. `pipeline.intext_citations --force` reuses existing TEI and
uses `processed.pdf` only when its receipt matches both byte hashes; it cannot
recover source-backed text from legacy TEI that lacks coordinates. Production
bundle and gold-subsample replays remain release acceptance work.

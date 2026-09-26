# Historical bibliography review — September 25–26

[Review receipt](review.json) records the audited `1.4.0.dev0` authority,
historical metadata and raw observations, current-code replay, formatter
results, and complete-refresh comparisons. Original authority and source
artifacts were preserved. These are historical-artifact replays, not results
from the fresh full-library candidate.

## #296: provenance and author precedence

The audited Mańko/Pugh metadata records BibTeX extraction, while its historical
authority row has no BibTeX import stamp. Replay restores the provenance and
the live formatter returns `provenance=bib`. Ordered Unicode authors were
already correct in this snapshot and remain correct after clean/incremental
replay. There is no stored reconciliation decision for this document; the
reported wrong-author merge mechanism is not reproduced or established.

## #314: Pugh split and acquisition ranking

The historical 59/23 split between the erroneous 1965 work and the canonical
1974 work becomes 21/61, with no overlap. All 38 moved observations contain
publication-year 1974 and the canonical title in retained raw citation text.
This review inspected those raw strings; it did not freshly inspect all 38
citing PDFs. All 90,180 raw observations are preserved without alteration.

The 21 remaining observations include OCR-damaged titles, omitted words,
mixed-script authors, and publication years after the journal rather than in
the author prefix. They are enumerated in the receipt. Twelve carry an explicit
possible-year-conflict warning; nine do not. The ghost is still ranked **13th**
in missing-work leads, with 21 citing papers and an aggregate warning. This
does not establish correct acquisition ranking. Pugh1990 and Alvarinoetal1990
remain unresolved using old extraction, despite their successful fresh-source
pilots. The fresh candidate must be reviewed separately.

Clean and incremental initial citation graphs agree, as do their top twelve
missing-work results. Agreement is not accuracy: the list still includes two
Kramp 1961 entries under different identities, among other review leads.

## First-refresh semantic equivalence fails

The initial replay had 88,818 citation edges; its first authority refresh
changed that graph. A later complete refresh of the retained 88,817-edge
database passed (job `27530129`), but that only demonstrates stability after
the initial transition. Its first reporting attempt (`27522972`) failed in
the harness because it mixed SQLite Row objects and tuples; that attempt's
database and logs are preserved separately.

Clean-cycle job `27537172` reproduced the first-transition discrepancy through
the complete authority/reconciliation sequence: **one edge disappears**, with
no added edges and all raw observations unchanged:

- Citing document: `188c66a35702`, work `10.2174/092986712803833308`.
- Cited work: `corpus:edwards|2000|edwards l hessinger d a toxicon 2000 38 `.
- The edge disappears during authority refresh; subsequent reconciliation
  does not restore it.

The cause and correct handling of that edge require investigation. This is
an unmet release equivalence gate, even though later unchanged refreshes are
stable. No graph rows were manually repaired to produce these results.

Raw scripts, databases and logs are retained under
[`scratch/v15-bouchet-20260925`](../../../scratch/v15-bouchet-20260925).

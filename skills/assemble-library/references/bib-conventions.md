# Keys, filenames, shelves — and the trap that has already corrupted a library

## Conventions

| Authors | Key | PDF filename |
|---|---|---|
| One | `Smith1998` | `Smith1998.pdf` |
| Two | `SmithJones1998` | `Smith_Jones1998.pdf` |
| Three or more | `Smithetal1998` | `Smithetal1998.pdf` |

Collisions take `a` / `b` / `c` suffixes: `Smith1998a`, `Smith1998b`.

PDFs are sharded into `library/<letter>/` by first-author surname initial, which
keeps a directory listing usable at a few thousand files. Files with no known
bibliographic record live in `library/orphans/` — kept, not discarded, because a
PDF whose metadata you have not yet worked out is still a paper you hold.

`file = {...}` appears **only** where a PDF is actually on disk. That field is
what a corpus build matches on, and it is what makes the bib a superset of the
library rather than a wish list indistinguishable from an inventory.

## The trap: positional suffixes silently rebind PDFs

Read this before writing key generation. It is not hypothetical — it happened,
and it currently blocks bib regeneration in an existing library.

**The `a`/`b`/`c` suffix is positional.** It falls out of an entry's rank among
its same-stem, same-year siblings. So merging a single duplicate shifts the
letter for every later sibling.

For entries with a real author that is survivable — you get a rename, you notice.
For author-less entries it is not. In the viburnum library roughly 165 entries
share the stem `Anon`, dozens of them with a PDF on disk named for one particular assignment,
and a shift silently rebinds `Anon2004.pdf` onto a *different paper*. Because
`file = {...}` is what the build matches on, wrong metadata then propagates into
every citation downstream — and nothing downstream can detect it, because the
entry looks fine and the PDF opens.

What makes this worse than an ordinary off-by-one: the merge that triggers it is
usually an *improvement*. Viburnum's `build_bib.py` refuses to write on 7 files
precisely because a title-cleaning fix merged 5 duplicate pairs that had looked
distinct. Fixing your data reshuffles your filenames.

### The guard, and why it is not the fix

That library's `build_bib.py` now **refuses to write** a bib in which any PDF
would attach to a different paper than the current bib records — naming the
affected files and leaving everything untouched. That is the right failure mode
and it saved the library. But it converts silent corruption into a hard block:
the committed bib is correct, and regeneration is what is now impossible.

**Do not inherit the bug in order to inherit the guard.** A new library should
be built so the guard has nothing to catch.

### Building it right

Derive stems from something **stable** rather than positional. Two approaches,
and the second is the honest one:

1. **A per-entry anchor** — DOI, else OpenAlex id, else an author/year/title
   signature. Note that none of these is invariant on its own: when a duplicate
   pair merges, the survivor inherits its twin's DOI *and* its author list, so an
   entry recorded as author-less with no DOI answers to none of its former
   identities.
2. **A committed lockfile mapping filename → work.** Explicit, reviewable, and
   survives any amount of metadata churn. Prefer this.

A third measure worth taking early: **backfill authors rather than tolerating an
`Anon` population**. Author-less entries are usually not anonymous — the harvest
record simply lacked an author list. Most of viburnum's were Donoghue Lab papers.
Filling those in retires the `Anon` stem entirely and removes the population where
the trap does its damage.

## Verifying a binding actually holds

Independent of key generation, a mis-resolved reference lands a real PDF of the
*wrong* paper under a plausible filename, and nothing downstream notices.

Compare each bib title against the extracted text of its PDF. Viburnum found
three such cases in 506 comparable papers by bag-of-words overlap against the
first ~2000 characters of `text.json` — one PDF was a *Kalina cycle* waste-heat
paper filed as a chloroplast genome. Restrict the comparison to Latin-script text
so translated titles do not all score zero.

This needs a built corpuscle, so it is a **post-build** check rather than part of
assembly — run it after the first `corpus run` and treat it as part of the first
improvement lap.

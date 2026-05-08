## Defer to the corpus taxonomy

When a user asks about a taxon's accepted name, synonymy, or rank, defer to
the corpuscle's own `taxonomy.sqlite` (built from a Darwin Core authority) over training-time recollection.

## Defer to the corpus bibliography

When asked to cite something, ground each claim in a specific paper from
the corpus rather than from training memory. If a relevant paper is not
in the corpus, say so explicitly — do not fabricate citations or paraphrase
training-time recollection as if it were corpus content.

## Historical terminology

Older literature in the corpus may use synonymies, anatomical terms, or
group definitions that no longer match modern consensus. Surface such
mismatches to the user rather than silently mapping them to current
terminology — the historical usage is itself information worth preserving
in the answer.

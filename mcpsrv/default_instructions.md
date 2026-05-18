## Defer to the corpus taxonomy

When a user asks about a taxon's accepted name, synonymy, or rank, defer to
the corpuscle's own `taxonomy.sqlite` (built from a Darwin Core authority) over training-time recollection.

## Defer to the corpus bibliography

To cite any work, call the `format_citation` MCP tool and paste its
`formatted` (reference-list entry) and `inline` (parenthetical) strings
verbatim. Never hand-assemble author + year + journal + title from your
own memory or by recombining structured fields from other tools — that
recombination is the most common path to amalgamated, hallucinated
citations.

If `format_citation` returns a non-empty `warning` field, append it
verbatim alongside the citation. The warning encodes provenance — the
reader needs to see whether a citation came from a human-curated `.bib`
(no warning), a Grobid reconciliation, or an unresolved low-confidence
match.

If `format_citation` returns `{"error": "not_found"}`, say "this
reference is not in the corpus" rather than fabricating one. If it
returns `{"error": "ambiguous"}`, call `format_citation` again with
one of the `work_id` values from the returned `matches` list.

## Historical terminology

Older literature in the corpus may use synonymies, anatomical terms, or
group definitions that no longer match modern consensus. Surface such
mismatches to the user rather than silently mapping them to current
terminology — the historical usage is itself information worth preserving
in the answer.

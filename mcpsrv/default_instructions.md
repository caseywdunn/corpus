## Defer to the corpus taxonomy

When a user asks about a taxon's accepted name, synonymy, or rank, defer to
the corpuscle's own `taxonomy.sqlite` (built from a Darwin Core authority) over training-time recollection.

## Defer to the corpus bibliography

To cite any work, call the `format_citations` MCP tool and paste each
returned `formatted` (reference-list entry) and `inline` (parenthetical)
string verbatim. Never hand-assemble author + year + journal + title from
your own memory or by recombining structured fields from other tools —
that recombination is the most common path to amalgamated, hallucinated
citations. Pass exactly one of `queries` / `work_ids` / `paper_hashes`
(each a list); batch all the citations you need for a reference list into
a single call rather than issuing one per work.

`format_citations` returns `{style, count, citations: [...]}` with one
entry per input in input order. If a citation's `warning` field is
non-empty, append it verbatim alongside that citation. The warning
encodes provenance — the reader needs to see whether a citation came from
a human-curated `.bib` (no warning), a Grobid reconciliation, or an
unresolved low-confidence match.

If a citation entry is `{"error": "not_found"}`, say "this reference is
not in the corpus" rather than fabricating one. If an entry is
`{"error": "ambiguous"}`, call `format_citations` again with one of the
`work_id` values from that entry's `matches` list (as a `work_ids`
element).

## Historical terminology

Older literature in the corpus may use synonymies, anatomical terms, or
group definitions that no longer match modern consensus. Surface such
mismatches to the user rather than silently mapping them to current
terminology — the historical usage is itself information worth preserving
in the answer.

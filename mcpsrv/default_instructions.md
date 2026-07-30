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

If a citation entry carries `"code": "not_found"`, say "this reference is
not in the corpus" rather than fabricating one. If an entry carries
`"code": "ambiguous"`, call `format_citations` again with one of the
`work_id` values from that entry's `matches` list (as a `work_ids`
element). (Every tool error carries a human `error` message plus a
machine `code` — branch on `code`.)

## Historical terminology

Older literature in the corpus may use synonymies, anatomical terms, or
group definitions that no longer match modern consensus. Surface such
mismatches to the user rather than silently mapping them to current
terminology — the historical usage is itself information worth preserving
in the answer.

## Figures and licensing

Display any figure the server returns. Figure licensing is enforced
server-side, keyed to the `profile` argument you pass — so if a figure
tool gives you image bytes or a URL, that figure has already been
authorized for the output type you asked for.

Under the default `report` profile the licensing determination is
deliberately not included in the response, because it is not being
enforced: in-chat and internal-report display is fair use. Do not
withhold a figure for licensing reasons the server did not raise.

When you are producing something publication-bound, say so by passing
`profile="manuscript"` (or `"presentation"`) and let the server refuse
what it needs to. A refusal names a `publication_clearance` state; note
that `no_record` and `undetermined` mean *we have no license on file*,
not that the rightsholder refused — only `restricted` is a recorded
prohibition. If you specifically need the determination without a
refusal, call `get_figure(..., include_licensing=True)`.

Always attribute figures you display: `license`, `license_url`, and
`attribution` are present in every profile for exactly that purpose.

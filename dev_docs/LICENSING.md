# Figure licensing + the publishable gate

How `corpus` reasons about whether a figure can be reused in a derived
publication, and what the `publishable` flag exposed by the MCP server
actually means. Backs [#51](https://github.com/caseywdunn/corpus/issues/51).

## Default policy

A figure is `publishable=true` when **either**:

1. The parent work carries an explicit `license` whose value is in the
   set treated as reusable: `public-domain`, `CC0-1.0`, `CC-BY-1.0`
   through `CC-BY-4.0`, `CC-BY-SA-3.0`, `CC-BY-SA-4.0`. Or
2. The parent work carries no license but its `year` is older than
   `licensing.pd_cutoff_years` years before the current year (default
   95, configurable per corpuscle).

A figure is `publishable=false` when:
- The license is `all-rights-reserved`, `publisher-permission`, or
  `unknown` (the explicit-not-publishable vocabulary), **or**
- No license is set and the year is too recent for the configured PD
  cutoff (or year is missing).

`publishable=null` happens when the license string is unrecognized
(neither in the publishable set nor in the not-publishable set).
Conservative readers should treat null as not-publishable.

## License sources

`works.license_source` records how the publishable decision was made:

| `license_source` | meaning |
|---|---|
| `bibtex` | the operator wrote `license = {...}` in the BibTeX file |
| `age_based_pd` | no license set; year falls outside the PD cutoff |
| `unknown` | no license, no year, or year too recent |

Both the raw fields (`license`, `license_url`) and the derived
`publishable` boolean are exposed via the MCP `get_figure` tool, so
clients with non-default jurisdictions or fair-use claims can re-derive.

## Jurisdiction caveat

The default `pd_cutoff_years: 95` matches the **US copyright rule**
(life-of-author + 70 years, with a hard 95-year cap on works
published before the rule changed). It does **not** cover all
jurisdictions:

- **EU**: life-of-author + 70 years; for older works without a known
  death year, 120 years is a safer default. Override per-corpuscle:
  ```yaml
  licensing:
    pd_cutoff_years: 120
  ```
- **Other jurisdictions**: consult local law or a librarian; the
  cutoff is a single integer per corpuscle so it's easy to tune.

This file is a starting point, not legal advice. Operators
publishing derived works are responsible for checking the rights
that apply in their jurisdiction.

## The publishable gate on `get_figure_image`

`mcpsrv.tools.figures.get_figure_image` raises a structured
`ValueError` when `publishable=false`:

> figure not publishable: license='unknown' (source='unknown'). The
> image is not returned to avoid downstream copyright issues. Read
> get_figure(<paper_hash>, <figure_id>) for the raw license fields if
> your jurisdiction / use case differs, or restart the server with
> --allow-unpublishable for local rights-holder cases.

The `--allow-unpublishable` server flag bypasses the gate
unconditionally. Use only when:

- You hold the rights to the figures (e.g. you authored the papers), **or**
- You're operating under fair use / fair dealing for a specific purpose
  and have made that determination yourself, **or**
- The corpus runs strictly locally and no derived publication will
  carry the figures.

Never set `--allow-unpublishable` on a public-facing deploy.

## Curating license fields

Use the BibTeX round-trip:

```bash
corpus bib export -o my_corpus.bib   # current works → BibTeX
$EDITOR my_corpus.bib                # add `license = {CC-BY-4.0}` etc.
corpus bib import my_corpus.bib      # apply edits back
```

License values follow [SPDX short identifiers](https://spdx.org/licenses/)
plus the small custom vocabulary above. `licenseurl` (flat) round-trips
to `works.license_url` (snake_case).

## Known limitations

- **Reprint chains** — a 1850 plate reproduced in a 2020 paper has the
  inner figure in the public domain even when the wrapper isn't. v0.3
  doesn't model this; the figure inherits the wrapper's license.
  Operators with reprint-heavy corpora should curate explicit license
  values at the figure level once the model supports per-figure
  overrides (out of scope for v0.3).
- **Photo-of-PD-work** — museum photographs of pre-1930 specimens may
  or may not carry separate copyright depending on jurisdiction. Not
  encoded; use explicit license values when known.
- **Per-figure overrides** — license fields apply to whole works.
  No per-figure `license` field exists yet; defer until a real case
  forces it.
- **Age computation** — uses the year field on the work, which can be
  Grobid-mis-parsed for older scans. Curate the year via BibTeX
  round-trip when accuracy matters for the cutoff decision.

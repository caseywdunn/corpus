# Figure licensing

How `corpus` reasons about whether a figure can be reused in a derived
publication, and what the server tells a client about it. Backs
[#51](https://github.com/caseywdunn/corpus/issues/51) and
[#154](https://github.com/caseywdunn/corpus/issues/154).

**The gate decides, not the client.** If a figure tool returns pixels or a
URL for the `profile=` you passed, that figure is cleared for that use —
display it. The clearance *determination* is only included in a response
when a strict profile is in force, or on an explicit
`get_figure(..., include_licensing=True)`.

## What a client sees: `publication_clearance`

Five states, reported as `publication_clearance` (#154 §2):

| state | meaning |
|---|---|
| `public_domain` | asserted public-domain, or past the PD cutoff |
| `licensed_open` | an explicit open license (CC-BY family) |
| `restricted` | terms are on file and forbid republication |
| `undetermined` | a license string was recorded but not recognized — usually a typo |
| `no_record` | nothing on file. **The common case, and not evidence of restriction.** |

Only `restricted` is a recorded prohibition. The distinction matters
because in a real corpus the overwhelming majority of works are
`no_record`: 86% of the served siphonophore corpus, with *not one* work
asserted `all-rights-reserved`. A single boolean collapsed those into the
same value, which read as a prohibition and caused clients to withhold
figures the server had cleared.

## The stored boolean

`works.publishable` remains the storage column and the input to the strict
gate — for *refusing*, "could not establish" and "explicitly restricted"
warrant the same conservative answer. It is no longer sent to clients.

A work is `publishable=1` when **either**:

1. The parent work carries an explicit `license` whose value is in the
   set treated as reusable: `public-domain`, `CC0-1.0`, `CC-BY-1.0`
   through `CC-BY-4.0`, `CC-BY-SA-3.0`, `CC-BY-SA-4.0`. Or
2. The parent work carries no license but its `year` is older than
   `licensing.pd_cutoff_years` years before the current year (default
   95, configurable per corpuscle).

It is `publishable=0` when:
- The license is `all-rights-reserved`, `publisher-permission`, or
  `unknown` (the explicit-not-publishable vocabulary), **or**
- No license is set and the year is too recent for the configured PD
  cutoff (or year is missing).

`publishable=NULL` happens when the license string is unrecognized —
neither vocabulary matched. That surfaces as `undetermined`, and the
authority build now logs a warning naming the offending string, because a
typo'd `license = {CC-BY 4.0}` (space, not hyphen) otherwise blocks
figures exactly like a real refusal, silently.

## License sources

`works.license_source` records how the publishable decision was made:

| `license_source` | meaning |
|---|---|
| `bibtex` | the operator wrote `license = {...}` in the BibTeX file |
| `age_based_pd` | no license set; year falls outside the PD cutoff |
| `unknown` | no license, no year, or year too recent |

`license`, `license_url`, and a server-computed `attribution` string go
out in **every** profile — captions need them. The clearance state does
not; see the note at the top.

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

## The figure-licensing gate, keyed to the output profile (#101)

The figure-licensing gate is governed by the active **output profile**,
which a client passes per call as `profile=` (see `mcpsrv/profiles.py`):

- **`report`** (the default) — *permissive*: `get_figure_image`,
  `get_figure_url`, and `get_figure_roi_image` return the figure
  regardless of clearance, on the basis that momentary in-chat display is
  fair use. No clearance determination is attached, precisely so it isn't
  mistaken for a permission check.
- **`manuscript`** / **`presentation`** — *strict*: the same three tools
  refuse, naming the state that caused it, e.g.

  > figure withheld under profile 'manuscript' — publication_clearance=
  > 'no_record': no license on file and the work is not old enough to be
  > age-based public domain; this is an ABSENCE of evidence, not a refusal
  > by the rightsholder. … For in-chat display request profile='report'.

All three pixel-returning tools enforce the gate; `get_figure_roi_image`
did not until #154 §3, which made it a way around the other two.
`tests/test_freeze_contract.py` now asserts every pixel-returning tool
takes a `profile`.

Selection is a **per-call client/session property**, not a server or
corpus setting — a shared SSE server serves chat, internal-report, and
manuscript sessions at once, so each call states its own profile. The
server's `--default-profile` only sets the fallback for calls that omit
`profile=` (default `report`); a startup warning fires when that
fallback is permissive. When `get_figure_url` issues a URL it encodes
the resolved profile into the URL, so the subsequent HTTP fetch enforces
the same policy.

> **Publication-bound work must pass `profile="manuscript"` per call.**
> The permissive default means a client that forgets will receive
> uncleared figures. The remedy is to state the profile, not to
> second-guess the server's answer.

(Replaces the v0.5 `--allow-unpublishable` server flag, which was a
single global toggle and could not distinguish concurrent sessions.)

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

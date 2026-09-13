# Retrieval ethics

These are constraints on how the library is built, not suggestions to weigh
against coverage. Carry their substance into the generated library's own
`CONTRIBUTING.md`, so the next person to extend the harvest inherits them.

## The rules

**Use documented, public APIs and the OA endpoints publishers themselves
advertise** — OpenAlex, Crossref, Unpaywall, Semantic Scholar, Europe PMC, BHL.

**Honour `robots.txt`, terms of service, HTTP 429 / `Retry-After`, and
documented rate limits.** Identify the client with a real contact address.

**Do not** attempt to defeat bot protection, rotate or spoof a User-Agent to
evade a block, use shadow-library mirrors (Sci-Hub and similar), or use anyone's
institutional credentials or proxy session.

**A 403 is a final answer, not a puzzle.** Record the failure and move on.

## Why this costs less than it looks

The honest handling is already the precedent in both existing libraries, and the
resulting libraries are good.

BHL's `partpdf` and `itempdf` endpoints return 403 to any non-interactive client
whatever you send, so viburnum's `fetch_pdfs.py` **deliberately does not call
them** — trying only wasted requests. Those entries land in the bib with a
`bhl_part` id, a `rights` statement and no `file`. Publisher 403s (MDPI, Wiley,
OUP, Allen Press) are the dominant fetch failure and are simply not engineered
around.

Expect roughly **a quarter to a third of attempts to succeed**. Every failure is
a want-list entry: the bib records the paper with its DOI or URL, and a human
fetches it by hand with institutional access. That is a working system, not a
degraded one — and it is why the bib-is-a-superset design exists.

One route nobody has taken: BHL's page-image API (`/pageimage/<PageID>`) does
serve, so assembling article PDFs page by page is possible within the rules. It
is not implemented anywhere, and it is a real option rather than a workaround.

## The contact address

Crossref, Unpaywall and OpenAlex all offer materially higher rate limits to
clients that identify themselves — the "polite pool". Use it.

**Read the address from the environment; never hardcode one.** Viburnum
hardcodes a `MAILTO` constant at the top of four scripts, and its
`CONTRIBUTING.md` has to tell anyone forking the library to go change all four —
because using someone else's address for a rate-limit courtesy is rude. Generated
templates should read it from an environment variable and fail with a clear
message when unset.

## API keys — the user's own, via environment only

Never write a key into the library repo. Check which are set, report what each
unlocks, and degrade gracefully to keyless endpoints.

| Variable | Unlocks |
|---|---|
| `BHL_API_KEY` | Biodiversity Heritage Library ([get one](https://www.biodiversitylibrary.org/getapikey.aspx)) |
| a contact-email variable | Crossref / Unpaywall / OpenAlex polite pools |
| others as the clade needs | e.g. Semantic Scholar |

Reporting which keys were absent is part of finishing the run: it tells the user
what a second pass would add, and distinguishes "this literature does not exist
in the open" from "you did not have the key for it".

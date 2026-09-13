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

**The variable is `CORPUS_CONTACT_EMAIL`.** Name it in the generated library's
readme and in every script that makes a request, so there is exactly one answer
to "what do I set?".

```bash
export CORPUS_CONTACT_EMAIL="you@example.edu"
```

Tell the user to set it **before the harvest starts**, and say what it buys: not
politeness in the abstract, but materially higher rate limits on three of the
services the harvest leans on hardest. A harvest run without it is slower and
more likely to hit a 429.

It also feeds the User-Agent, which is the other half of identifying yourself:

```python
CONTACT = os.environ.get("CORPUS_CONTACT_EMAIL")
if not CONTACT:
    raise SystemExit(
        "CORPUS_CONTACT_EMAIL is not set. Public indexes give identified "
        "clients much higher rate limits, and an unidentified harvest is both "
        "slower and ruder.\n"
        '  export CORPUS_CONTACT_EMAIL="you@example.edu"'
    )
UA = f"corpus-library/1.0 (mailto:{CONTACT})"
```

**Fail with that message rather than falling back to a default.** There is no
safe default: a placeholder address is a lie to the service, and someone else's
real address sends them the consequences of your traffic.

**Never hardcode it.** Viburnum hardcodes a `MAILTO` constant at the top of four
scripts, and its `CONTRIBUTING.md` has to tell anyone forking the library to go
change all four — because using someone else's address for a rate-limit courtesy
is rude. Reading one environment variable in one place is the fix; four
constants and a note in the docs is the bug.

It is a contact address, not a credential: it belongs in the generated readme as
an instruction, and never in a committed file as a value.

## API keys — the user's own, via environment only

| Variable | Unlocks | Secret? |
|---|---|---|
| `BHL_API_KEY` | Biodiversity Heritage Library — the pre-1930 literature | **yes** |
| `CORPUS_CONTACT_EMAIL` | Crossref / Unpaywall / OpenAlex polite pools | no — see above |
| others as the clade needs | e.g. Semantic Scholar | usually |

### Obtaining a BHL key

Request one at <https://www.biodiversitylibrary.org/getapikey.aspx>. It is free,
tied to an email address, and arrives by email. BHL's API v3 takes it as an
`apikey=` query parameter.

Confirm the current process on that page rather than trusting this paragraph —
it sits behind bot protection, which this skill will not work around, so it
cannot be checked automatically.

**Without it, the pre-1930 half of the library does not happen.** BHL is the only
practical route to article-level records for that material, so for a clade with
a deep historical literature this is the difference between a corpus that starts
in 1990 and one that starts in 1758. Say that to the user rather than treating
the key as optional polish.

### Setting one

Any of these work; none of them put the key in a file that could be committed:

```bash
# A dedicated secrets file, sourced when you need it — easiest to revoke
mkdir -p ~/.config/corpus
read -rs BHL_API_KEY   # typed, not echoed, not in shell history
printf 'export BHL_API_KEY=%s\n' "$BHL_API_KEY" > ~/.config/corpus/secrets.env
chmod 600 ~/.config/corpus/secrets.env

source ~/.config/corpus/secrets.env && python scripts/harvest_bhl.py
```

Or a `.env` at the library root — already gitignored in this project's
convention, and what corpus itself uses for `ANTHROPIC_API_KEY` — or an export
in a shell profile if the key is long-lived and the machine is yours alone.

### The rules for handling one

**Never write a key into the library repo**, including in a script default, a
committed `.env`, a config file, or a comment. Check that it is set; do not
store it.

**Never paste a key into a chat or an issue.** It lands in the transcript, and
transcripts get logged, shared and pasted elsewhere. When an agent needs a key,
the answer is a file it can read, never a message it can quote.

**Degrade gracefully rather than failing.** A missing key should skip that
source with a named warning, not abort a harvest that can still do most of its
work. `CORPUS_CONTACT_EMAIL` is the exception — it is not a secret, costs
nothing to set, and every service the harvest touches rewards it.

**Report which keys were absent when the run finishes.** It tells the user what
a second pass would add, and it distinguishes "this literature does not exist in
the open" from "you did not have the key for it" — two very different answers
that an unreported skip makes indistinguishable.

**If a key leaks, revoke it rather than rotating the file.** Request a new one at
the URL above; a key in a transcript stays in that transcript.

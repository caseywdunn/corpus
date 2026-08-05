# API stability policy (1.0)

What "frozen" means for the MCP surface, what counts as a breaking
change, and how something gets removed.

Deferred out of v0.6 by design: the surface was still moving, and a
policy written over a moving surface documents intentions rather than
commitments. It is written now because 1.0 is the version other people
build against.

The short version: **the 38 tools, their parameter names, and the shape
of what they return are a contract. Adding to them is a minor release;
changing or removing anything in them requires a major.**

## What is covered

The contract is the MCP surface a client sees over the wire:

| Covered | Specifically |
| --- | --- |
| **Tool set** | The 38 names in `EXPECTED_TOOLS` (`tests/test_freeze_contract.py`) |
| **Parameter names + meanings** | Including which are optional and what their defaults do |
| **Response field names** | A field present in a 1.0 response stays present, with the same meaning |
| **Pagination convention** | `limit` is the single-knob name; only the seven documented multi-cap dossier/graph tools use `max_*` |
| **Error payload** | Every error carries a human `error` and a machine `code` drawn from `ERROR_CODES` |
| **Output profiles** | `report` / `manuscript` / `presentation` and what each permits |
| **Licensing enforcement** | Every tool that returns figure pixels accepts `profile` and consults the shared gate |

Four of those seven rows are enforced by `tests/test_freeze_contract.py`
rather than by review, which is the only reason to believe them. The test
is the contract; this document explains it.

## What is not covered

Deliberately outside the guarantee, because freezing them would freeze
the pipeline's ability to improve:

- **On-disk artifacts.** `metadata.json`, `figures.json`, `text.json` and
  the rest carry their own `ARTIFACT_SCHEMA_VERSION` (`pipeline/__init__.py`),
  independent of the API version. A corpuscle is not a public API; it is
  an intermediate the server reads.
- **Extracted *content*.** Better OCR, better chunking, and better taxon
  matching change what the tools return, and that is the point of the
  project. v1.0's re-OCR change rewrote stored text for most scanned
  papers. **A stable API does not promise stable data.** Clients that
  need reproducibility should pin a built corpuscle, not a corpus
  version.
- **The `corpus` CLI.** Covered by ordinary semver on the distribution,
  but it is an operator tool, not the integration surface. Flags may be
  added freely and renamed at a minor with a deprecation warning.
- **Anything under `mcpsrv/` that is not a registered tool.** Helpers,
  index internals, and `_`-prefixed functions are private, whatever their
  import path suggests.
- **`dev_docs/`.** Internal.

## Breaking vs. additive

**Additive — ships in a minor release (1.1, 1.2, …):**

- A new tool.
- A new *optional* parameter whose default preserves existing behavior.
- A new field in a response object.
- A new `code` value in the error vocabulary. Clients must treat an
  unrecognized `code` as a generic failure rather than crashing — the
  vocabulary is closed at any given version but is expected to grow.
- A new output profile.

**Breaking — requires a major release (2.0):**

- Removing or renaming a tool, a parameter, or a response field.
- Changing a default, since a client that omitted the parameter gets
  different behavior without changing a line.
- Narrowing a response: dropping a field, or making a previously
  unconditional field conditional.
- Tightening a type or an accepted value range.
- Changing what an error `code` means, as distinct from adding one.
- Changing a tool's semantics while keeping its signature — the worst
  kind, because nothing fails loudly.

The last two are the ones that get argued about. The test is not "did the
signature change" but **"can a correct 1.0 client get a worse answer, or
no answer, without changing its code?"** If yes, it is breaking.

### The worked example

v1.0 removed `publishable` from every figure response and replaced it
with `publication_clearance`, which reports five states where
`publishable` collapsed "the rightsholder forbade this" and "we have no
license record" into a single `0` (#154). That is squarely breaking under
the rules above.

It shipped in 1.0 precisely because 1.0 is the last moment it is cheap.
The same change proposed against 1.4 waits for 2.0, or ships as an
additive field alongside the old one with `publishable` deprecated. That
is the whole practical content of the freeze.

## Deprecation path

Nothing in the covered surface disappears without going through this:

1. **Announce.** The tool or field keeps working. Its docstring — which
   is the tool description every MCP client reads, so this is the place
   clients actually see it — gains a `Deprecated since X.Y:` line naming
   the replacement. A CHANGELOG entry under **Deprecated** says the same.
2. **Wait.** At least one full minor release, and at least 90 days. A
   deprecation announced in 1.3 cannot be removed in 1.4 if 1.4 lands a
   month later.
3. **Remove at the next major**, listed in the CHANGELOG under
   **Changed (breaking, MCP surface)** with the migration in the same
   entry.

Server-side, a deprecated field stays populated. Emitting `null` from a
field a client still reads is a removal wearing a deprecation's clothes.

There is no runtime deprecation *warning* channel: MCP has no out-of-band
place to put one, and stuffing it into the response would itself change
the response shape. The docstring and the CHANGELOG are the channel.

## How a client should pin

`bundle_info` returns `server_version` (the code) and `bundle_version`
(the corpuscle), and they move independently — the same 1.2.0 server can
serve a corpuscle built by 1.0.1. A client that cares should read both
and check `server_version` against the major it was written for.

Pin the corpuscle, not just the server, when reproducibility matters. The
server guarantees the *shape* of an answer; only a fixed corpuscle
guarantees the answer.

## Changing the contract

Editing `EXPECTED_TOOLS` is editing the public API. The test says so and
will fail until the set is updated deliberately. The process:

1. Decide whether it is additive or breaking, by the rules above.
2. If breaking, it does not go on `dev` — it goes on the branch for the
   next major, or it is redesigned as an additive change.
3. Update `EXPECTED_TOOLS` and any other assertion in
   `tests/test_freeze_contract.py` in the *same commit* as the code.
4. Add the CHANGELOG entry, with a migration note if breaking.
5. Regenerate the table in [MCP_TOOLS.md](MCP_TOOLS.md).

A change that reaches review with `EXPECTED_TOOLS` edited and no
CHANGELOG entry is incomplete, not nearly done.

## Related

- [MCP_TOOLS.md](MCP_TOOLS.md) — the surface itself
- [../CONTRIBUTING.md](../CONTRIBUTING.md) — release ritual and versioning
- [../CHANGELOG.md](../CHANGELOG.md) — what actually changed
- `tests/test_freeze_contract.py` — the enforced half of this document

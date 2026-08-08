"""CITATION.cff must stay schema-valid, or Zenodo stops archiving releases.

This file was invalid from the day it was added (2026-05-10) until v1.1.1,
and the consequence was invisible: Zenodo's webhook receiver answers
``202`` on receipt and only *then* builds the deposition from
CITATION.cff, so an invalid file fails after GitHub has already recorded
a successful delivery. Seven releases — v0.3.0 through v1.1.0 — were
never archived, and the concept DOI kept resolving to v0.2.0, so every
citation of this software pointed at a three-month-old version.

Two violations did it, and both are cheap to catch:

    year: 2026                    # CFF 1.2.0 has no top-level `year`
    license: "(see repository)"   # must be an SPDX identifier

Scope: this is a targeted subset of CFF 1.2.0, not a full schema
validator — it checks the shape of failure that actually occurred plus
the required-field floor. `cffconvert --validate -i CITATION.cff` is the
real thing; run it if you are making structural changes. A subset test
that runs on every push is worth more here than a complete one nobody
installs, because the failure mode is silence rather than an error.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).parent.parent
CFF_PATHS = [REPO_ROOT / "CITATION.cff", REPO_ROOT / "pipeline" / "CITATION.cff"]

# CFF 1.2.0 top-level keys. A key outside this set is what broke Zenodo.
_ALLOWED_TOP_LEVEL = frozenset({
    "abstract", "authors", "cff-version", "commit", "contact",
    "date-released", "doi", "identifiers", "keywords", "license",
    "license-url", "message", "preferred-citation", "references",
    "repository", "repository-artifact", "repository-code", "title",
    "type", "url", "version",
})

_REQUIRED = ("cff-version", "message", "title", "authors")

# Not the full SPDX list — the point is to reject free prose like
# "(see repository)", which is what was there. An identifier is short,
# has no spaces, and is drawn from a constrained alphabet.
_SPDX_CHARS = set(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.-+"
)


def _load(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


@pytest.mark.parametrize("path", CFF_PATHS, ids=lambda p: str(p.name))
def test_cff_parses_as_yaml(path):
    data = _load(path)
    assert isinstance(data, dict), f"{path} did not parse to a mapping"


@pytest.mark.parametrize("path", CFF_PATHS, ids=lambda p: str(p.name))
def test_cff_has_no_unknown_top_level_keys(path):
    """The exact failure: `year: 2026` is not a CFF key (`date-released` is)."""
    unknown = sorted(set(_load(path)) - _ALLOWED_TOP_LEVEL)
    assert not unknown, (
        f"{path.name}: {unknown} not valid CFF 1.2.0 top-level key(s). "
        f"An unknown key fails schema validation, and Zenodo answers 202 "
        f"before it discovers that — so the release looks fine and is "
        f"never archived. Did you mean `date-released` or `version`?"
    )


@pytest.mark.parametrize("path", CFF_PATHS, ids=lambda p: str(p.name))
def test_cff_has_required_keys(path):
    data = _load(path)
    missing = [k for k in _REQUIRED if k not in data]
    assert not missing, f"{path.name}: missing required CFF key(s) {missing}"


@pytest.mark.parametrize("path", CFF_PATHS, ids=lambda p: str(p.name))
def test_cff_license_looks_like_an_spdx_identifier(path):
    """The other failure: `license: "(see repository)"` is free prose."""
    data = _load(path)
    if "license" not in data:
        pytest.skip("no license key (optional in CFF)")
    licenses = data["license"]
    for lic in [licenses] if isinstance(licenses, str) else licenses:
        assert lic and not (set(lic) - _SPDX_CHARS), (
            f"{path.name}: license {lic!r} is not an SPDX identifier. "
            f"CFF constrains this to the SPDX enum — free text fails "
            f"validation. Use e.g. MIT, Apache-2.0, GPL-3.0-or-later."
        )


def test_cff_copies_are_identical():
    """#173 — `pipeline/CITATION.cff` is a hand-maintained copy with nothing
    enforcing sync. Nothing still enforces *why* it exists, but at least a
    silent divergence now fails here."""
    a, b = (p.read_bytes() for p in CFF_PATHS)
    assert a == b, (
        "CITATION.cff and pipeline/CITATION.cff have diverged. They are "
        "kept byte-identical by hand (#173); update both."
    )

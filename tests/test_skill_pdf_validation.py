"""A paywall page must never land in a library looking like a paper (#178).

`skills/assemble-library/scripts/fetch_pdfs.py` guards the one boundary in the
harvest where a wrong answer is invisible downstream. Publisher paywalls return
HTTP 200 with an HTML body; if one is written to `library/` under a plausible
filename, a corpus build ingests it as the paper its bib entry names, extracts
nothing useful, and no later stage has any way to notice. The bib says it is
Smith 1998 and the file opens.

So the validator is tested here rather than left to the harvest. It is a pure
function over bytes precisely so this can run offline, with no network and no
downloads — which is also why these cases are fabricated rather than fetched.

What the boundary actually is, measured rather than assumed: `pdfinfo` catches
structural damage (truncation, a missing trailer, a broken xref) and is lenient
about semantic damage — a valid PDF of the *wrong paper* passes. That gap is
closed after a build by comparing bib titles against extracted text, and is
documented in the skill's `references/bib-conventions.md`. Do not add a case
here asserting otherwise; it would be asserting a guarantee this layer does not
make.
"""
from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "skills" / "assemble-library" / "scripts" / "fetch_pdfs.py"

pytestmark = pytest.mark.skipif(
    not SCRIPT.exists(), reason="assemble-library skill not present"
)


def _load():
    """Import the script by path — skills ship as templates, not as a package."""
    spec = importlib.util.spec_from_file_location("_fetch_pdfs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _minimal_pdf(pad: int = 30_000) -> bytes:
    """A structurally valid single-page PDF, padded past the size floor."""
    body = (
        b"%PDF-1.4\n"
        b"1 0 obj<</Type/Catalog/Pages 2 0 R>>endobj\n"
        b"2 0 obj<</Type/Pages/Kids[3 0 R]/Count 1>>endobj\n"
        b"3 0 obj<</Type/Page/Parent 2 0 R/MediaBox[0 0 612 792]>>endobj\n"
    )
    xref_at = len(body)
    body += b"xref\n0 4\n0000000000 65535 f \n"
    for off in (9, 52, 108):
        body += b"%010d 00000 n \n" % off
    body += b"trailer<</Size 4/Root 1 0 R>>\nstartxref\n"
    body += str(xref_at).encode() + b"\n%%EOF\n"
    return body + b"%" + b"p" * pad + b"\n"


@pytest.fixture(scope="module")
def validate():
    return _load().is_acceptable_pdf


@pytest.fixture(scope="module")
def floor():
    return _load().MIN_BYTES


@pytest.mark.parametrize(
    "name, blob_factory",
    [
        ("html interstitial", lambda n: b"<!DOCTYPE html><html>Access denied</html>" + b" " * n),
        ("html after whitespace", lambda n: b"\n\n  <html>paywall</html>" + b" " * n),
        ("binary junk", lambda n: b"\x00\x01not a pdf" + b"\x00" * n),
    ],
)
def test_non_pdf_bodies_are_rejected(validate, floor, name, blob_factory):
    """The realistic paywall shapes — all HTTP 200, none of them papers."""
    ok, reason = validate(blob_factory(floor))
    assert not ok, f"{name} was accepted"
    assert reason


def test_html_is_named_as_such(validate, floor):
    """The reason matters: 'html interstitial' tells an operator it is a
    paywall, where 'not a PDF' reads like a broken link and invites a retry."""
    ok, reason = validate(b"<!DOCTYPE html><html>Sign in</html>" + b" " * floor)
    assert not ok
    assert "html" in reason.lower()


def test_undersized_is_rejected_even_with_a_pdf_header(validate):
    """A few kB of 'access denied' can still start with %PDF-."""
    ok, reason = validate(b"%PDF-1.4 tiny")
    assert not ok
    assert "too small" in reason


def test_a_real_pdf_is_accepted(validate, tmp_path):
    blob = _minimal_pdf()
    path = tmp_path / "ok.pdf"
    path.write_bytes(blob)
    ok, reason = validate(blob, path)
    assert ok, f"valid PDF rejected: {reason}"


@pytest.mark.skipif(
    subprocess.run(["which", "pdfinfo"], capture_output=True).returncode != 0,
    reason="pdfinfo not installed",
)
@pytest.mark.parametrize(
    "name, damage",
    [
        ("truncated body", lambda b: b[:60] + b"p" * 30_000),
        ("trailer removed", lambda b: b.split(b"xref")[0] + b"p" * 30_000),
        ("broken xref root", lambda b: b.replace(b"/Root 1 0 R", b"/Root 9 0 R")),
    ],
)
def test_structurally_damaged_pdfs_are_rejected(validate, tmp_path, name, damage):
    """The truncated-download case, which is the common real failure."""
    blob = damage(_minimal_pdf())
    path = tmp_path / "bad.pdf"
    path.write_bytes(blob)
    ok, reason = validate(blob, path)
    assert not ok, f"{name} was accepted"
    assert "pdfinfo" in reason

"""Tests for reference extraction and metadata quality."""

import pytest


def test_min_reference_count(paper_references):
    """Pipeline extracts at least the expected number of references."""
    gt, refs_data = paper_references
    min_count = gt["references"].get("min_count")
    if min_count is None:
        pytest.skip("No min_count specified")
    actual = refs_data.get("total_references", len(refs_data.get("references", [])))
    assert actual >= min_count, f"Too few references: {actual} < {min_count}"


def test_expected_references(paper_references):
    """Specific expected references are found."""
    gt, refs_data = paper_references
    expected = gt["references"].get("expected", [])
    if not expected:
        pytest.skip("No specific references to check")
    refs = refs_data.get("references", [])

    for exp in expected:
        def matches_ref(ref):
            if "year" in exp and ref.get("year") != exp["year"]:
                return False
            if "title_contains" in exp:
                title = (ref.get("title") or "").lower()
                if exp["title_contains"].lower() not in title:
                    return False
            if "authors_contain" in exp:
                ref_authors = " ".join(ref.get("authors") or []).lower()
                for a in exp["authors_contain"]:
                    if a.lower() not in ref_authors:
                        return False
            return True

        found = [r for r in refs if matches_ref(r)]
        desc = ", ".join(f"{k}={v}" for k, v in exp.items())
        assert found, f"Expected reference not found: {desc}"


def test_metadata_title(paper_metadata):
    """Paper title is correctly extracted."""
    gt, meta = paper_metadata
    expected = gt["metadata"].get("title")
    if expected is None:
        pytest.skip("No title in ground truth")
    actual = meta.get("title", "")
    assert expected.lower().strip() in actual.lower().strip(), (
        f"Title mismatch: expected {expected!r}, got {actual!r}"
    )


def test_metadata_year(paper_metadata):
    """Publication year is correctly extracted."""
    gt, meta = paper_metadata
    expected = gt["metadata"].get("year")
    if expected is None:
        pytest.skip("No year in ground truth")
    assert meta.get("year") == expected, (
        f"Year mismatch: expected {expected}, got {meta.get('year')}"
    )


def test_metadata_authors(paper_metadata):
    """Expected authors appear in the extracted author list."""
    gt, meta = paper_metadata
    expected_authors = gt["metadata"].get("authors_contain", [])
    if not expected_authors:
        pytest.skip("No authors in ground truth")
    actual_surnames = [a.get("surname", "").lower() for a in meta.get("authors", [])]
    for exp in expected_authors:
        surname = exp.get("surname", "").lower()
        assert surname in actual_surnames, (
            f"Author {surname!r} not found in {actual_surnames}"
        )

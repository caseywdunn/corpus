"""Tests for text extraction quality."""

import pytest


def test_must_contain(paper_text):
    """Extracted text contains expected domain phrases."""
    gt, text_data = paper_text
    text = text_data.get("text", "").lower()
    for phrase in gt["text"].get("must_contain", []):
        assert phrase.lower() in text, f"Expected phrase not found: {phrase!r}"


def test_must_not_contain(paper_text):
    """Extracted text does not contain known garbage strings."""
    gt, text_data = paper_text
    text = text_data.get("text", "").lower()
    for phrase in gt["text"].get("must_not_contain", []):
        assert phrase.lower() not in text, f"Garbage phrase found: {phrase!r}"


def test_min_length(paper_text):
    """Extracted text meets minimum length threshold."""
    gt, text_data = paper_text
    min_len = gt["text"].get("min_length")
    if min_len is None:
        pytest.skip("No min_length specified")
    actual = len(text_data.get("text", ""))
    assert actual >= min_len, f"Text too short: {actual} < {min_len}"

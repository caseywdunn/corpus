"""Caption-binding evidence exposed by the thin MCP compatibility layer."""

from mcpsrv.tools.figures import _caption_evidence_fields


def test_explicit_build_evidence_is_preserved():
    fields = _caption_evidence_fields({
        "caption_text": "Figure 1. Description.",
        "page": 2,
        "caption_page": 3,
        "caption_status": "uncertain",
        "caption_confidence": "low",
        "caption_page_distance": 1,
        "caption_kind": "prose_caption",
    })

    assert fields == {
        "caption_status": "uncertain",
        "caption_confidence": "low",
        "caption_page_distance": 1,
        "caption_kind": "prose_caption",
    }


def test_legacy_structural_same_page_caption_is_bound():
    fields = _caption_evidence_fields({
        "caption_text": "Figure 1. Description.",
        "caption_source": "docling_caption_link",
        "page": 2,
        "caption_page": 2,
    })

    assert fields["caption_status"] == "bound"
    assert fields["caption_confidence"] == "high"
    assert fields["caption_page_distance"] == 0
    assert fields["caption_kind"] == "unknown"


def test_legacy_cross_page_caption_is_not_presented_as_ordinary():
    fields = _caption_evidence_fields({
        "caption_text": "Plate IV.",
        "caption_source": "heuristic_proximity",
        "page": 4,
        "caption_page": 5,
    })

    assert fields["caption_status"] == "uncertain"
    assert fields["caption_confidence"] == "low"
    assert fields["caption_page_distance"] == 1
    assert fields["caption_kind"] == "unknown"


def test_legacy_blank_caption_is_unbound():
    fields = _caption_evidence_fields({
        "caption_text": "",
        "page": 4,
    })

    assert fields["caption_status"] == "unbound"
    assert fields["caption_confidence"] is None
    assert fields["caption_page_distance"] is None
    assert fields["caption_kind"] is None

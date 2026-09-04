"""Caption-selection regressions from the gold corpuscle (#195).

These are deliberately geometry-level fixtures. They preserve the Docling
shape that caused a real figure on page 12 to inherit ``FIGURE 9`` from page
13: the correct ``FIGURE 8`` label was the second provenance span of a text
item whose first span was running prose on page 11.
"""
from types import SimpleNamespace as NS

from pipeline.figures import extract_caption_info


def _bbox(left, bottom, right, top):
    return NS(l=left, b=bottom, r=right, t=top)


def _prov(page, coords, charspan=None):
    return NS(page_no=page, bbox=_bbox(*coords), charspan=charspan)


def _text(text, page, coords, *, label="text", charspan=None, prov=None):
    spans = prov if prov is not None else [_prov(page, coords, charspan)]
    return NS(text=text, label=label, prov=spans)


def _picture(page, coords, *, captions=None):
    return NS(prov=[_prov(page, coords, (0, 0))], captions=captions or [])


class _Ref:
    def __init__(self, target):
        self.target = target

    def resolve(self, _document):
        return self.target


def test_page_spanning_text_exposes_the_correct_figure_label():
    """The Hosia 2024 page-12 failure: Fig. 8 must not become Fig. 9."""
    body = "Palpons continue from the preceding page."
    merged = f"{body} FIGURE 8"
    split_label = NS(
        text=merged,
        label="text",
        prov=[
            _prov(11, (305, 77, 546, 123), (0, len(body))),
            _prov(12, (61, 234, 85, 239), (len(body) + 1, len(merged))),
        ],
    )
    caption8 = _text(
        "Ontogenetic change in Nanomia cara nectophore shape.",
        12,
        (61, 216, 524, 232),
        label="caption",
    )
    label9 = _text("FIGURE 9", 13, (61, 506, 85, 511), label="section_header")
    caption9 = _text(
        "Composite illustration of mature Nanomia cara nectophore and bracts.",
        13,
        (61, 479, 521, 503),
        label="caption",
    )
    figure8 = _picture(12, (84, 247, 511, 754))
    figure9 = _picture(13, (48, 468, 546, 766))
    document = NS(
        texts=[split_label, caption8, label9, caption9],
        pictures=[figure8, figure9],
    )

    got = extract_caption_info(figure8, document)

    assert got["caption_text"].startswith("FIGURE 8.")
    assert "Ontogenetic change" in got["caption_text"]
    assert got["caption_page"] == 12
    assert got["caption_page_distance"] == 0
    assert got["caption_status"] == "bound"
    assert got["caption_confidence"] == "medium"
    assert got["caption_kind"] == "prose_caption"
    rejected = [c for c in got["caption_candidates"] if not c["chosen"]]
    assert any(
        c["caption_text"].startswith("FIGURE 9")
        and c["rejection_reason"] == "substantial_picture_on_caption_page"
        for c in rejected
    )


def test_adjacent_label_is_joined_to_a_structurally_linked_caption_body():
    label = _text("FIGURE 3.", 4, (60, 300, 90, 308), label="section_header")
    body = _text(
        "A complete descriptive caption for the extracted figure.",
        4,
        (60, 270, 500, 298),
        label="caption",
    )
    picture = _picture(4, (60, 320, 500, 700), captions=[_Ref(body)])
    document = NS(texts=[label, body], pictures=[picture])

    got = extract_caption_info(picture, document)

    assert got["caption_text"].startswith("FIGURE 3.")
    assert "complete descriptive caption" in got["caption_text"]
    assert got["caption_source"] == "docling_caption_link"
    assert got["caption_confidence"] == "high"
    assert "FIGURE 3.." not in got["caption_text"]


def test_structural_label_selects_its_slice_of_an_enumerating_caption():
    label = _text("Fig. 11.", 7, (146, 189, 180, 196), label="caption")
    grouped = _text(
        "Fig. 10. Colony and sections. Fig. 11. Eudoxid. "
        "Fig. 12. Isolated gonophore.",
        7,
        (100, 140, 296, 183),
        label="caption",
    )
    picture = _picture(7, (137, 205, 199, 371), captions=[_Ref(label)])
    document = NS(texts=[label, grouped], pictures=[picture])

    got = extract_caption_info(picture, document)

    assert got["caption_text"] == "Fig. 11. Eudoxid."
    assert got["caption_kind"] == "prose_caption"


def test_candidate_diagnostic_is_bounded_without_truncating_selected_caption():
    body_text = "A " + ("long descriptive caption " * 40)
    label = _text("FIGURE 2", 3, (60, 300, 90, 308))
    body = _text(body_text, 3, (60, 270, 500, 298), label="caption")
    picture = _picture(3, (60, 320, 500, 700))
    document = NS(texts=[label, body], pictures=[picture])

    got = extract_caption_info(picture, document)

    assert got["caption_text"].endswith("caption")
    assert len(got["caption_text"]) > 600
    assert len(got["caption_candidates"][0]["caption_text"]) == 600


def test_facing_page_candidate_survives_when_no_local_picture_competes():
    picture = _picture(4, (40, 60, 550, 760))
    label = _text("Plate IV.", 5, (50, 700, 110, 715))
    document = NS(texts=[label], pictures=[picture])

    got = extract_caption_info(picture, document)

    assert got["caption_text"] == "Plate IV."
    assert got["caption_page"] == 5
    assert got["caption_page_distance"] == 1
    assert got["caption_status"] == "uncertain"
    assert got["caption_confidence"] == "low"
    assert got["caption_kind"] == "bare_label"


def test_no_candidate_is_explicitly_unbound():
    picture = _picture(2, (40, 60, 550, 760))
    document = NS(texts=[_text("Running prose.", 2, (50, 20, 500, 40))], pictures=[picture])

    got = extract_caption_info(picture, document)

    assert got["caption_text"] == ""
    assert got["caption_status"] == "unbound"
    assert got["caption_confidence"] is None
    assert got["caption_kind"] is None
    assert got["caption_candidates"] == []

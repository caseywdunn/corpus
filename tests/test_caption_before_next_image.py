"""Issue-derived #336 layout; no source extraction or model run is claimed."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from PIL import Image

from pipeline.figure_passes import _pass25_annotate_figures, _pass3b_annotate_rois
from pipeline.figures import (
    classify_figure, dedupe_figures, expand_plate_figures, extract_caption_info,
    plate_legend_entries,
)
from tests.test_caption_association import _picture, _Ref, _text


def layout(*, linked_next=False, extra_caption=False):
    source = json.loads((Path(__file__).parent / "fixtures/figure_integrity/"
                         "caption_before_next_image.json").read_text())
    texts = [_text(t["text"], t["page"], t["bbox"], label="caption")
             for t in source["texts"]]
    pictures = [_picture(p["page"], p["bbox"], captions=(
        [_Ref(texts[p["caption_index"]])] if p["caption_index"] is not None else []
    )) for p in source["pictures"]]
    if linked_next:
        pictures[1].captions = [_Ref(texts[1])]
    document = SimpleNamespace(texts=texts, pictures=pictures)
    items = [{"docling_idx": p["docling_idx"], "page": p["page"], "bbox": p["bbox"],
              **extract_caption_info(picture, document)}
             for p, picture in zip(source["pictures"], pictures)]
    captions = [{"text": t["text"], "bbox": t["bbox"]} for t in source["texts"]]
    if extra_caption:
        captions.append({"text": "Fig. 5: Another separate plot.", "bbox": [50, 10, 550, 40]})
    return items, captions


@pytest.mark.parametrize("extra_caption", [False, True])
def test_separate_captions_do_not_clone_previous_image(extra_caption):
    # Adding a third caption must not turn conflicting ownership into proof.
    items, captions = layout(extra_caption=extra_caption)
    result = expand_plate_figures(items, {178: plate_legend_entries(captions)})
    assert result == items
    assert [item["figure_number"] for item in result] == ["3", None]
    assert result[1]["caption_status"] == "unbound"
    rejected = [c for c in result[0]["caption_candidates"]
                if c["rejection_reason"] == "separate_caption_with_following_picture"]
    assert {c["figure_number"] for c in rejected} == ({"4", "5"} if extra_caption else {"4"})
    assert all(not c["chosen"] and c["competing_docling_indices"] == [96] for c in rejected)
    before = json.dumps(result, sort_keys=True)
    assert json.dumps(expand_plate_figures(items, {178: plate_legend_entries(captions)}),
                      sort_keys=True) == before


def test_real_structural_link_can_bind_next_page_without_a_false_clone():
    items, captions = layout(linked_next=True)
    result = expand_plate_figures(items, {178: plate_legend_entries(captions)})
    assert len(result) == 2
    assert result[1]["figure_number"] == "4"
    assert result[1]["caption_status"] == "bound"
    assert result[1]["caption_page_distance"] == -1
    assert result[1]["caption_source"] == "docling_caption_link"
    assert not any("shares_image_with" in item for item in result)


@pytest.mark.parametrize("positive", ["plate_heading", "grouped_caption", "bare_host"])
def test_independent_plate_evidence_survives_an_uncaptioned_following_image(positive):
    items, captions = layout()
    if positive == "plate_heading":
        captions.insert(0, {"text": "PLATE XII", "bbox": [50, 720, 550, 750]})
    elif positive == "grouped_caption":
        captions = [{"text": "Fig. 3. Colony. Fig. 4. Eudoxid.", "bbox": [50, 60, 550, 110]}]
    else:
        items[0]["caption_text"] = "Fig. 3."
    result = expand_plate_figures(items, {178: plate_legend_entries(captions)})
    assert result[-1]["figure_number"] == "4"
    assert result[-1]["shares_image_with"] == 95


def test_next_page_furniture_is_not_a_competing_figure():
    items, captions = layout()
    items[1]["bbox"] = [0, 0, 5, 5]
    result = expand_plate_figures(items, {178: plate_legend_entries(captions)})
    assert result[-1]["shares_image_with"] == 95


def test_vision_receives_only_actual_panel_targets_not_false_figure4(tmp_path):
    items, captions = layout()
    items = expand_plate_figures(items, {178: plate_legend_entries(captions)})
    for item in items:
        item["figure_type"] = classify_figure(item)
    items = dedupe_figures(items)
    for item in items:
        item["figure_id"] = f"docling_{item['docling_idx']}"
        path = tmp_path / f"{item['figure_id']}.png"
        Image.new("RGB", (200, 200), "white").save(path)
        item.update(file_path=str(path), filename=path.name)
    figures_file, text_file = tmp_path / "figures.json", tmp_path / "text.json"
    figures_file.write_text(json.dumps({"figures": items}))
    text_file.write_text(json.dumps({"text": "The phylogeny is shown in Fig. 4."}))
    _pass25_annotate_figures(text_file, figures_file)
    calls = []

    class Backend:
        name = "vision:test"

        def detect_figure_panels(self, path, caption, labels):
            calls.append((path.name, caption, labels))
            return []

    _pass3b_annotate_rois(figures_file, Backend())
    assert len(calls) == 1
    assert calls[0][0] == "docling_95.png"
    assert calls[0][2] == ["A", "B", "C"]
    assert "Phylogeny" not in calls[0][1]
    records = json.loads(figures_file.read_text())["figures"]
    assert len(records) == 2
    assert next(record for record in records if record["figure_id"] == "docling_96")[
        "caption_status"] == "unbound"
    assert not any(record.get("plate_figures_from_caption") for record in records)

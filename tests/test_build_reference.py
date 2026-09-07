"""Release-reference comparisons must expose equal-count evidence corruption."""
import copy
import json

import pytest

from tools.qc.build_reference import ARTIFACTS, compare, main, snapshot


@pytest.fixture
def build(tmp_path):
    hd = tmp_path / "documents" / "a"
    hd.mkdir(parents=True)
    for name in ARTIFACTS:
        (hd / name).write_text(json.dumps({"schema_version": 1, "pipeline_version": "1.3.dev0"}))
    (hd / "references.json").write_text(json.dumps({"references": [
        {"xml_id": "b0", "title": "Active and passive factors affecting aggregations",
         "authors": ["K Kosobokova", "R Hopcroft"], "raw": "Source citation"}],
        "schema_version": 1}))
    (hd / "summary.json").write_text(json.dumps({"processing_summary": {
        "quality_flags": [{"gate": "zero_references_unexpected", "severity": "warn"}]}}))
    return tmp_path


def test_same_count_reference_corruption_requires_review(build):
    old = snapshot(build)
    path = build / "documents/a/references.json"
    refs = json.loads(path.read_text())
    refs["references"][0]["title"] = "Active and passive factors"
    refs["references"][0]["authors"] = ["K Hopcroft", "R"]
    path.write_text(json.dumps(refs))
    diff = compare(old, snapshot(build))
    assert diff["requires_review"]
    assert diff["changed"]["a"]["artifacts"] == ["references.json"]
    fields = diff["changed"]["a"]["references"][0]["fields"]
    assert fields["authors"]["before"] == ["K Kosobokova", "R Hopcroft"]
    assert fields["title"]["after"] == "Active and passive factors"


def test_unchanged_warning_is_reported_not_automatically_failed(build):
    report = snapshot(build)
    diff = compare(report, report)
    assert not diff["requires_review"]
    assert diff["unchanged_warning_documents"] == ["a"]


def test_unchanged_hard_failure_still_fails(build):
    report = snapshot(build)
    report["documents"]["a"]["quality_flags"][0]["severity"] = "error"
    assert compare(report, report)["requires_review"]


def test_additions_removals_missing_artifacts_and_bad_schema(build):
    old = snapshot(build)
    (build / "documents/a/text.json").unlink()
    new = snapshot(build)
    assert new["problems"] == ["a/text.json: missing"]
    assert compare(old, new)["requires_review"]
    moved = copy.deepcopy(old)
    moved["documents"]["b"] = moved["documents"].pop("a")
    diff = compare(old, moved)
    assert diff["added"] == ["b"] and diff["removed"] == ["a"]
    moved["schema_version"] = 999
    with pytest.raises(ValueError, match="resnapshot"):
        compare(old, moved)


def test_root_relocation_and_producer_stamp_do_not_hide_content(build, tmp_path_factory):
    import shutil
    path = build / "documents/a/figures.json"
    path.write_text(json.dumps({"pipeline_version": "old", "figures": [
        {"image_path": str(build / "documents/a/figures/image.png"), "caption_text": "caption"}]}))
    old = snapshot(build)
    other = tmp_path_factory.mktemp("relocated") / "build"
    shutil.copytree(build, other)
    path = other / "documents/a/figures.json"
    path.write_text(path.read_text().replace(str(build), str(other)).replace('"old"', '"new"'))
    assert not compare(old, snapshot(other))["requires_review"]
    path.write_text(path.read_text().replace('"caption"', '"wrong"'))
    assert compare(old, snapshot(other))["requires_review"]


def test_cli_never_overwrites_reference(build):
    path = build / "reference.json"
    assert main(["snapshot", str(build), "--out", str(path)]) == 0
    before = path.read_bytes()
    with pytest.raises(SystemExit) as exc:
        main(["snapshot", str(build), "--out", str(path)])
    assert exc.value.code == 2
    assert path.read_bytes() == before
    assert main(["compare", str(path), str(path)]) == 0


def test_binary_changes_are_distinguished_from_semantic_diffs(build):
    path = build / "documents/a/grobid.tei.xml"
    path.write_text("first generated ID")
    old = snapshot(build)
    path.write_text("another generated ID")
    diff = compare(old, snapshot(build))
    assert diff["binary_changes"] == {"a": ["grobid.tei.xml"]}
    assert not diff["requires_review"]


def test_same_count_index_and_pixel_changes_require_review(build):
    old = snapshot(build)
    new = copy.deepcopy(old)
    new["indexes"]["vectors"] = {"changed": "same number of rows"}
    assert compare(old, new)["index_changes"]
    assert compare(old, new)["requires_review"]
    new = copy.deepcopy(old)
    new["documents"]["a"]["figure_pixels"] = {"figure.png": {"rgba_sha256": "changed"}}
    assert compare(old, new)["changed"]["a"]["figure_pixels"]
    assert compare(old, new)["requires_review"]

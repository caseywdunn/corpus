"""Frozen retrieval workflows distinguish source hits, assistance and missing review."""
import copy
import json
from pathlib import Path

import pytest

from tools.qc.retrieval import (
    capture_local, compare_reports, digest, enrich_rows, evaluate, expected_calls,
    sample_independent, score_rows, target_match, validate_manifest,
)

MANIFEST = Path(__file__).parents[1] / "dev_docs/examples/siphonophore_retrieval_evaluation.json"


def target(paper="abc", page=4):
    return {"id": "diagnosis", "paper_hash": paper, "pages": [page],
            "all_text": ["distinctive diagnostic character"], "grade": 2,
            "review_status": "source_verified", "source": {"physical_page": page, "sha256": paper}}


def manifest():
    return {"queries": [{"id": "q", "group": "independent",
                          "call": {"query": "taxon diagnosis", "k": 5, "with_text": True},
                          "targets": [target()]}],
            "acceptance": {"groups": {"independent": {"hit5": 1, "hit10": 1}}}}


def capture(m, rows):
    return {"manifest_sha256": digest(m), "identity": {"label": "synthetic test; no ranking claim"},
            "runs": [{"query_id": "q", "mode": "unassisted", "variant": variant,
                      "call": call, "rows": rows} for variant, call in expected_calls(m["queries"][0]).items()]}


def hit():
    return {"paper_hash": "abc", "text": "A DISTINCTIVE diagnostic character.", "source_pages": [4],
            "table_refs": [], "table_metadata_available": True}


def test_source_anchors_require_content_and_reject_wrong_document_or_known_page():
    assert target_match(target(), hit())
    assert target_match(target(), dict(hit(), source_pages=[]))  # explicitly counted legacy route
    for row in (dict(hit(), paper_hash="other"), dict(hit(), source_pages=[8]),
                dict(hit(), text="Taxon diagnosis ........ 4")):
        assert not target_match(target(), row)
    assert not target_match(dict(target(), review_status="pending"), hit())


def test_exact_seven_audit_queries_and_three_controls_remain_frozen():
    m = json.loads(MANIFEST.read_text())
    validate_manifest(m)
    actual = [q["call"] for q in m["queries"] if q["group"] == "audit"]
    expected = [
        ("How can I distinguish Nanomia bijuga from Nanomia cara?", 5, None),
        ("Nanomia bijuga versus Nanomia cara diagnostic characters", 10, None),
        ("Wie unterscheiden sich Nanomia cara und Nanomia bijuga?", 5, None),
        ("How many species of Physalia are recognized and what distinguishes them?", 5, None),
        ("key to the species of Erenna distinguishing characters nectophore tentilla", 5, "da370f1e0434"),
        ("Apolemia lanosa diagnosis", 5, "c9e7e8ae50a2"),
        ("Key to Erenna species based on their Type A bracts", 5, "da370f1e0434"),
    ]
    assert actual == [dict(query=q, k=k, with_text=True, **({"paper_hash": p} if p else {})) for q, k, p in expected]
    assert [q["call"]["query"] for q in m["queries"] if q["group"] == "prose_control"] == [
        "pneumatophore structure and gas gland function", "Erenna bioluminescent lures fish prey capture",
        "somatocyst function in calycophoran siphonophores"]
    assert len([q for q in m["queries"] if q["group"] == "historical_control"]) == 2
    assert "not the audited" in m["baseline_warning"]


def test_hit_at_five_and_ten_use_distinct_ranks_and_do_not_assert_precision():
    m = manifest()
    rows = [dict(hit(), text="unjudged comparative passage")] * 5 + [hit()]
    report = evaluate(m, capture(m, rows))
    assert report["status"] == "fail"
    assert report["gates"][0]["hit_rate"] == 0
    assert report["gates"][1]["hit_rate"] == 1
    result = report["queries"][-1]
    assert result["unjudged_rows"] == 5
    assert result["known_target_coverage"] == 1


def test_missing_source_review_or_runs_block_release_acceptance():
    m = manifest()
    m["queries"][0]["targets"] = []
    assert evaluate(m, capture(m, [hit()]))["status"] == "blocked"
    m = manifest()
    c = capture(m, [hit()])
    c["runs"].pop(0)
    report = evaluate(m, c)
    assert report["status"] == "blocked"
    assert report["missing_unassisted_runs"] == [{"query_id": "q", "variant": "original"}]
    c = capture(m, [{"error": "embedding unavailable"}])
    assert evaluate(m, c)["status"] == "blocked"
    c = capture(m, [hit()])
    c["runs"][0]["rows"] = [{"error": "original call failed"}]
    report = evaluate(m, c)
    assert all(g["status"] == "pass" for g in report["gates"])
    assert report["status"] == "blocked"
    assert report["failed_unassisted_runs"] == [{"query_id": "q", "variant": "original"}]


def test_assisted_recovery_never_relabels_unassisted_failure():
    m = manifest()
    c = capture(m, [])
    c["runs"].append({"query_id": "q", "mode": "assisted", "variant": "paper_scoped",
                       "call": {"query": "taxon diagnosis", "paper_hash": "abc", "k": 5}, "rows": [hit()]})
    report = evaluate(m, c)
    assert report["status"] == "fail"
    assert report["queries"][-1]["hit"] is True
    assert all(g["hits"] == 0 for g in report["gates"])
    c["runs"][-1]["mode"] = "unassisted"
    with pytest.raises(ValueError, match="differ"):
        evaluate(m, c)


def test_query_mutation_manifest_drift_and_duplicate_runs_are_rejected():
    m = manifest()
    c = capture(m, [])
    c["runs"][0]["call"]["paper_hash"] = "abc"
    with pytest.raises(ValueError, match="differ"):
        evaluate(m, c)
    c = capture(m, [])
    m["queries"][0]["targets"][0]["all_text"] = ["changed after seeing rank"]
    with pytest.raises(ValueError, match="different manifest"):
        evaluate(m, c)
    c = capture(m, [])
    c["runs"].append(copy.deepcopy(c["runs"][0]))
    with pytest.raises(ValueError, match="duplicate"):
        evaluate(m, c)


def test_diversity_and_repeated_table_metrics_distinguish_unknown_provenance():
    q = manifest()["queries"][0]
    table = dict(hit(), table_refs=["#/tables/1"])
    rows = [table, table, dict(table, text="different logical table row"), dict(hit(), paper_hash="xyz")]
    result = score_rows(q, rows, 5)
    assert result["unique_documents"] == 2 and result["dominant_document_fraction"] == .75
    assert result["repeated_table_rate"] == .5
    assert result["duplicate_table_text_rate"] == .25
    result = score_rows(q, [dict(r, table_metadata_available=False) for r in rows], 5)
    assert result["repeated_table_rate"] is None
    assert result["repeated_table_rate_known_lower_bound"] == .5


def test_candidate_gates_and_independent_improvement_are_separate():
    m = manifest()
    passing = evaluate(m, capture(m, [hit()]))
    failing = evaluate(m, capture(m, []))
    assert compare_reports(failing, passing)["independent_improvement_demonstrated"]
    unchanged = compare_reports(passing, passing)
    assert unchanged["candidate_release_gates"] == "pass"
    assert not unchanged["independent_improvement_demonstrated"]
    changed = dict(passing, manifest_sha256="different")
    with pytest.raises(ValueError, match="same frozen"):
        compare_reports(passing, changed)


def test_independent_source_sampling_is_reproducible_deduplicated_and_ungraded(tmp_path):
    m = manifest()
    m["queries"][0]["group"] = "audit"
    for paper in ("abc", "bcd", "cde", "def", "efg"):
        directory = tmp_path / "documents" / paper
        directory.mkdir(parents=True)
        chunks = []
        for page in (4, 5, 6):
            for split in (0, 1):
                chunks.append({"chunk_id": f"{page}_{split}", "section_type": "diagnosis",
                               "treatment_context": {"status": "resolved", "name": f"Species {paper}"},
                               "source_items": [{"item_ref": f"#/texts/{page}", "page": page}]})
        (directory / "chunks.json").write_text(json.dumps({"chunks": chunks}))
    result = sample_independent(m, tmp_path, 10, 320)
    assert result == sample_independent(m, tmp_path, 10, 320)
    sample = result["independent_sampling"]
    assert sample["eligible_units"] == 14  # duplicate split chunks and audit page excluded
    assert sample["selected"] == 10 and sample["selected_papers"] == 5
    assert all(not (c["paper_hash"] == "abc" and 4 in c["pages"]) for c in sample["selection"])
    queries = [q for q in result["queries"] if q["group"] == "independent"]
    assert len(queries) == 10 and all(q["targets"] == [] for q in queries)
    assert all("paper_hash" not in q["call"] for q in queries)  # unassisted corpus-wide retrieval
    assert m["queries"][0]["targets"]  # input manifest not mutated
    queries[0]["targets"] = [target(queries[0]["source_review_candidate"]["paper_hash"],
                                    queries[0]["source_review_candidate"]["pages"][0])]
    repeated = sample_independent(result, tmp_path, 10, 320)
    assert repeated["independent_sampling"] == sample


def test_enrichment_uses_materialized_source_metadata_without_inference(tmp_path):
    directory = tmp_path / "documents" / "abc"
    directory.mkdir(parents=True)
    artifact = {"treatment_context_policy": "source_treatments_v2", "chunks": [{
        "chunk_id": "c1", "source_items": [{"page": 1, "source_page": 4}],
        "tables": [{"table_ref": "#/tables/2"}],
    }]}
    (directory / "chunks.json").write_text(json.dumps(artifact))
    rows = enrich_rows([dict(hit(), chunk_id="c1")], tmp_path, {})
    assert rows[0]["source_pages"] == [4]
    assert rows[0]["table_refs"] == ["#/tables/2"]
    assert rows[0]["table_metadata_available"]
    legacy = enrich_rows([dict(hit(), paper_hash="legacy", chunk_id="c1")], tmp_path, {})
    assert legacy[0]["table_metadata_available"] is False
    artifact["chunks"][0] = {"chunk_id": "c1", "source_items": [], "treatment_context": {"status": "unknown"}}
    (directory / "chunks.json").write_text(json.dumps(artifact))
    naive = enrich_rows([dict(hit(), chunk_id="c1")], tmp_path, {})
    assert naive[0]["table_metadata_available"] is False


def test_independent_source_exclusion_uses_original_physical_pages(tmp_path):
    m = manifest()
    m["queries"][0]["group"] = "audit"
    directory = tmp_path / "documents" / "abc"
    directory.mkdir(parents=True)
    row = {"chunk_id": "c1", "section_type": "diagnosis",
           "treatment_context": {"status": "resolved", "name": "Species example"},
           "source_items": [{"page": 1, "source_page": 4}]}
    (directory / "chunks.json").write_text(json.dumps({"chunks": [row]}))
    result = sample_independent(m, tmp_path, 10, 320)
    assert result["independent_sampling"]["eligible_units"] == 0
    row["source_items"][0]["source_page"] = 8
    (directory / "chunks.json").write_text(json.dumps({"chunks": [row]}))
    result = sample_independent(m, tmp_path, 10, 320)
    assert result["independent_sampling"]["selection"][0]["pages"] == [8]


def test_committed_verified_anchors_exist_in_the_referenced_source_fragments():
    root = Path(__file__).parents[1]
    m = json.loads(MANIFEST.read_text())
    for query in m["queries"]:
        for label in query["targets"]:
            source = label["source"]
            fixture = json.loads((root / source["fixture"]).read_text())
            if "case" in source:
                region = fixture[source["case"]]
                sha = region["sha256"]
                texts = [i.get("text", "") for i in region["items"]]
            else:
                sha = fixture["source"]["sha256"]
                texts = [i.get("text", "") for i in fixture["document"].get("texts", [])]
                texts += [c["text"] for t in fixture["document"].get("tables", []) for c in t["data"]["table_cells"]]
            assert source["sha256"] == sha and label["paper_hash"] == sha[:12]
            assert target_match(label, {"paper_hash": sha[:12], "text": "\n".join(texts)})


def test_local_capture_uses_exact_tool_calls_and_records_identity(tmp_path, monkeypatch):
    from mcpsrv import app, indexes
    from mcpsrv.tools import chunks
    calls = []

    class Index:
        def __init__(self, output):
            self.papers = {"abc": {}}
            self.bundle_manifest = {"bundle_version": "1.2.1"}
            self._embedding_identity = {"model": "test-only"}

        def load(self):
            return 1

    def query(**kwargs):
        calls.append(kwargs)
        return [hit()]

    monkeypatch.setattr(indexes, "CorpusIndex", Index)
    monkeypatch.setattr(chunks, "get_chunks_for_topic", query)
    original_index = app._INDEX
    m = manifest()
    result = capture_local(m, tmp_path, "older-reference", "reference")
    assert calls == list(expected_calls(m["queries"][0]).values())
    assert result["identity"]["bundle_manifest"]["bundle_version"] == "1.2.1"
    assert result["identity"]["historical_audit_equivalence"].startswith("not_asserted")
    assert result["manifest_sha256"] == digest(m)
    assert app._INDEX is original_index

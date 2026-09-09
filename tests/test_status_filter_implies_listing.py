"""`--filter-*` on its own lists the affected papers (#169), and the
naive-chunker fallback is visible in the rollup (#168).

Both are the CLI lying by omission. The status report printed
`List affected papers with: corpus status --filter-gate <name>`, and
running exactly that reprinted the whole report unfiltered, because the
flag only took effect alongside `--list-hashes`. And a chunking fallback
that collapses retrieval quality logged once per paper and appeared in no
summary at all.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline import status as status_mod
from pipeline.stages import _run_quality_gates
from pipeline.status import _GATE_INFO, aggregate, render_text


def _paper(documents_dir: Path, h: str, flags=(), failures=()):
    hd = documents_dir / h
    hd.mkdir(parents=True, exist_ok=True)
    (hd / "summary.json").write_text(json.dumps({
        "pdf_hash": h,
        "processing_summary": {
            "stage_timings": [],
            "stage_failures": list(failures),
            "quality_flags": list(flags),
        },
    }))


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    docs = tmp_path / "documents"
    docs.mkdir()
    _paper(docs, "aaa")
    _paper(docs, "bbb", flags=[{"gate": "empty_text", "severity": "error"}])
    _paper(docs, "ccc", flags=[{"gate": "empty_text", "severity": "error"}])
    _paper(docs, "ddd", flags=[{"gate": "low_text_density", "severity": "warn"}])
    return tmp_path


# --- #169 ---------------------------------------------------------------


def test_a_filter_alone_lists_the_matching_papers(output_dir, capsys, monkeypatch):
    """The exact invocation the report's hint promises."""
    monkeypatch.setattr(
        "sys.argv",
        ["corpus-status", str(output_dir), "--filter-gate", "empty_text"],
    )
    assert status_mod.main() == 0
    printed = capsys.readouterr().out.split()
    assert sorted(printed) == ["bbb", "ccc"]


def test_the_hint_names_an_invocation_that_works(output_dir):
    """The report used to tell the reader to run a flag that did nothing
    on its own. Whatever the hint says now has to be the working form."""
    rendered = render_text(aggregate(output_dir / "documents"))
    hint = next(line for line in rendered.splitlines()
                if "List affected papers with" in line)
    assert "--filter-gate" in hint
    assert "--list-hashes" not in hint


def test_list_hashes_without_a_filter_is_unchanged(output_dir, capsys, monkeypatch):
    monkeypatch.setattr("sys.argv",
                        ["corpus-status", str(output_dir), "--list-hashes"])
    assert status_mod.main() == 0
    assert sorted(capsys.readouterr().out.split()) == ["bbb", "ccc", "ddd"]


def test_a_filter_that_cannot_apply_says_so_instead_of_vanishing(
    output_dir, monkeypatch, caplog,
):
    """`--json` does not filter. Silently ignoring the flag is the defect
    this issue is about, so the other direction has to be loud too."""
    monkeypatch.setattr(
        "sys.argv",
        ["corpus-status", str(output_dir), "--json",
         "--filter-gate", "empty_text"],
    )
    with caplog.at_level("WARNING"):
        assert status_mod.main() == 0
    assert "--filter-gate" in caplog.text
    assert "--json" in caplog.text
    assert "ignored" in caplog.text


# --- #168 ---------------------------------------------------------------


def _chunked_by(tmp_path: Path, chunker: str, n_chunks: int = 1) -> Path:
    hd = tmp_path / "doc"
    hd.mkdir(parents=True, exist_ok=True)
    (hd / "text.json").write_text(json.dumps(
        {"text": "body text " * 200, "pages": 2}))
    (hd / "chunks.json").write_text(json.dumps({
        "chunker": chunker,
        "total_chunks": n_chunks,
        "chunks": [{"chunk_id": f"chunk_{i}", "text": "x" * 50}
                   for i in range(n_chunks)],
    }))
    return hd


def test_the_naive_chunker_fallback_raises_a_quality_flag(tmp_path):
    """A 2-page paper chunked to 1 window instead of 16 — the run exits 0
    and every other gate passes."""
    gates = _run_quality_gates(_chunked_by(tmp_path, "naive_char_window"))
    hit = [g for g in gates if g["gate"] == "naive_chunker_fallback"]
    assert len(hit) == 1
    assert hit[0]["severity"] == "error"
    assert "naive_char_window" in hit[0]["detail"]


def test_the_hybrid_chunker_raises_nothing(tmp_path):
    gates = _run_quality_gates(_chunked_by(tmp_path, "hybrid_chunker", 16))
    assert [g for g in gates if g["gate"] == "naive_chunker_fallback"] == []


def test_an_artifact_with_no_chunker_recorded_raises_nothing(tmp_path):
    """Older builds predate the field; absence is not evidence of a
    fallback."""
    hd = _chunked_by(tmp_path, "hybrid_chunker", 16)
    data = json.loads((hd / "chunks.json").read_text())
    del data["chunker"]
    (hd / "chunks.json").write_text(json.dumps(data))
    gates = _run_quality_gates(hd)
    assert [g for g in gates if g["gate"] == "naive_chunker_fallback"] == []


def test_the_gate_is_explained_in_the_status_report():
    """An unexplained gate name in the rollup sends an operator to the
    source. Every other gate carries its remedy; so must this one."""
    severity, description = _GATE_INFO["naive_chunker_fallback"]
    assert severity == "error"
    assert "corpus prefetch" in description


def test_the_fallback_is_countable_across_a_corpus(tmp_path):
    """The original cause degraded *every* paper on a host following the
    offline recipe, so what matters is that the rollup can count it."""
    docs = tmp_path / "documents"
    docs.mkdir()
    for h in ("aaa", "bbb"):
        _paper(docs, h, flags=[{"gate": "naive_chunker_fallback",
                                "severity": "error"}])
    _paper(docs, "ccc")
    rollup = aggregate(docs)
    assert rollup["quality_flags"]["naive_chunker_fallback"] == 2
    assert "naive_chunker_fallback" in render_text(rollup)

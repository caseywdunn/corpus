"""Why a re-run is doing work, rolled up by reason (#80).

Per-stage implicit resume is correct and opaque. The drift computation
itself landed with v1.3's fingerprint work — `configuration_drift` and
`source_input_drift` already derive, from the same receipts resume uses,
which stages would re-run and why. What was missing was a way to read it:
the renderer printed one line per affected document, each repeating the
same handful of reasons. On the 699-document Viburnum corpuscle that is
699 lines; on the 1,775-document siphonophore one, 1,775.

#281 was that gap at its worst — the GPU vision phase silently
re-extracting every document, turning a 1.5-hour phase into a projected
35 — and finding it took hours of log archaeology for an answer that fits
on one line: `docling_extraction` re-runs on all of them because
`pipeline_version` changed.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from pipeline.status import render_drift_rollup


def test_a_reason_affecting_every_document_says_so_in_one_line():
    """The case that was hardest to see and matters most."""
    differences = {
        f"doc{i:04d}": {"docling_extraction": ["pipeline_version"]}
        for i in range(699)
    }
    out = render_drift_rollup(differences, total=699)
    assert out.splitlines() == [
        "  docling_extraction: pipeline_version — all of 699 documents"
    ]


def test_reasons_are_ordered_by_how_many_documents_they_affect():
    differences = {
        "aaa": {"docling_extraction": ["pipeline_version"],
                "text_chunking": ["config"]},
        "bbb": {"docling_extraction": ["pipeline_version"]},
        "ccc": {"docling_extraction": ["pipeline_version"]},
    }
    lines = [l for l in render_drift_rollup(differences, 3).splitlines()
             if not l.startswith("      ")]
    assert lines[0].startswith("  docling_extraction: pipeline_version")
    assert lines[1].startswith("  text_chunking: config")


def test_a_partial_reason_names_the_documents():
    """When only some documents drift, which ones is the useful detail —
    and it is exactly the case where a rollup alone would hide the answer."""
    differences = {
        "aaa": {"pdf_preparation": ["ocrlang"]},
        "bbb": {"pdf_preparation": ["ocrlang"]},
    }
    out = render_drift_rollup(differences, total=100)
    assert "  pdf_preparation: ocrlang — 2 of 100 documents" in out
    assert "      aaa, bbb" in out


def test_the_document_list_is_capped_and_says_how_many_it_hid():
    differences = {f"doc{i}": {"scan_detection": ["config"]} for i in range(12)}
    out = render_drift_rollup(differences, total=100, detail_limit=5)
    assert "scan_detection: config — 12 of 100 documents" in out
    assert "+7 more" in out
    listed = [l for l in out.splitlines() if l.startswith("      ")][0]
    names = [t for t in listed.strip().split(", ") if not t.startswith("+")]
    assert len(names) == 5, names


def test_a_full_scope_reason_does_not_list_documents():
    """Listing all 699 hashes under "all of 699" is the thing being
    replaced."""
    differences = {f"doc{i}": {"quality_gates": ["config"]} for i in range(699)}
    out = render_drift_rollup(differences, total=699)
    assert not any(l.startswith("      ") for l in out.splitlines())


def test_one_document_reads_as_singular():
    out = render_drift_rollup({"aaa": {"scan_detection": ["config"]}}, total=50)
    assert "1 of 50 document" in out
    assert "documents" not in out


def test_no_drift_says_nothing_will_re_run():
    """Silence would read as "the check did not run"."""
    assert "no stage would re-run" in render_drift_rollup({}, total=699)


def test_several_reasons_for_one_stage_are_counted_separately():
    """A stage re-running for two different reasons is two facts, and
    fixing one may not stop the re-run."""
    differences = {"aaa": {"pdf_preparation": ["ocrlang", "keeppages"]}}
    out = render_drift_rollup(differences, total=1)
    assert "pdf_preparation: ocrlang" in out
    assert "pdf_preparation: keeppages" in out


# --- the integration point the issue names ------------------------------


def test_dry_run_explains_the_plan_it_prints():
    """`--dry-run` already showed *what* would run; the issue's own
    suggestion was to add *why*. Only on dry-run: the check takes ~7 s on
    699 documents, which is cheap for a plan and not free enough to put in
    front of every build."""
    import inspect

    from pipeline import cli
    src = inspect.getsource(cli._cmd_run)
    assert "if args.dry_run:" in src
    assert "_explain_resume(cfg, config_path)" in src


def test_explaining_a_plan_never_blocks_the_plan():
    """Read-only diagnosis. An unreadable receipt or an unresolvable config
    must be reported and stepped over, not turned into a failure."""
    import inspect

    from pipeline import cli
    src = inspect.getsource(cli._explain_resume)
    assert "except (OSError, ValueError, TypeError, RuntimeError)" in src
    assert "return" in src
    assert "sys.exit" not in src


def test_a_first_build_is_not_reported_as_drift(tmp_path, capsys):
    """Every stage runs on a corpuscle with no documents/ tree, and that is
    a first build rather than drift. Calling it drift would train operators
    to ignore the warning."""
    from pipeline import cli
    from pipeline.config_schema import CorpuscleConfig

    cfg = CorpuscleConfig(input_pdfs=tmp_path / "pdfs",
                          output_dir=tmp_path / "out")
    cli._explain_resume(cfg, tmp_path / "config.yaml")
    assert capsys.readouterr().out == ""

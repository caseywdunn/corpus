"""Grobid consolidation flags are reachable from config (#PLAN v1.2 §3).

`consolidateCitations` had never run in this project's history:
`grobid_client.process_fulltext` defaulted it to 0, `metadata.py` called that
method overriding nothing, and the `grobid:` config block exposed only `url`
and `disable`. So the setting existed and could not be changed.

The default stays off, but it is now a *measured* default. Same PDFs with the
flag off and on, against the gold corpuscle:

    Ahuja 2026            86.1% -> 88.9% DOIs    3.1s -> 6.4s
    Stepanjants 2014       0.0% ->  6.9%         2.0s -> 2.7s
    Bernstein 1934         0.0% ->  3.6%         2.3s -> 3.4s
    Benasso 1976           0.0% ->  0.0%         1.7s -> 2.6s

Six DOIs across 194 references, for 1.4x to 2x the Grobid time. The split is
by era, not by anything the flag controls: modern reference lists already
carry DOIs and CrossRef holds the rest, while the historical works this corpus
is mostly made of are not in it at all. Which is exactly why the flag is
exposed rather than hard-coded — the arithmetic belongs to the library.
"""
from __future__ import annotations

import inspect

import pytest

from pipeline.config_schema import CorpuscleConfig, validate_config


def test_defaults_are_header_on_citations_off():
    c = CorpuscleConfig()
    assert c.grobid.consolidate_header == 1
    assert c.grobid.consolidate_citations == 0


def test_citation_consolidation_can_be_turned_on():
    """The point of the change: it was unreachable before."""
    cfg = validate_config({"grobid": {"consolidate_citations": 1}})
    assert cfg.grobid.consolidate_citations == 1


def test_out_of_range_levels_are_rejected():
    with pytest.raises(Exception):
        validate_config({"grobid": {"consolidate_citations": 9}})
    with pytest.raises(Exception):
        validate_config({"grobid": {"consolidate_header": -1}})


def test_metadata_passes_the_configured_levels_to_grobid():
    """Guard the wiring, which is where this was broken. `process_fulltext`
    always had the parameters; nothing ever supplied them."""
    from pipeline import metadata
    src = inspect.getsource(metadata)
    assert "consolidate_citations=int(" in src
    assert "consolidate_header=int(" in src


def test_the_measured_default_is_recorded_next_to_it():
    """A default chosen by measurement is only useful if the measurement is
    findable from the default."""
    from pipeline import grobid_client
    doc = grobid_client.GrobidClient.process_fulltext.__doc__
    assert "measured" in doc
    assert "194 references" in doc


def test_the_template_documents_both_flags():
    from pathlib import Path
    tpl = (Path(__file__).parent.parent / "pipeline" / "config.template.yaml").read_text()
    assert "consolidate_header" in tpl and "consolidate_citations" in tpl

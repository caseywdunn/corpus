"""Every exposed collection validates the same budget before corpus work."""
import inspect
import typing

import pytest

from mcpsrv import app
from mcpsrv.app import (
    MAX_COLLECTION_CHARS, MAX_COLLECTION_ITEMS, MAX_COLLECTION_ITEM_CHARS,
    _validate_collection,
)
from tests.test_freeze_contract import _fn, _registered


def collection_parameters():
    found = []
    for name, tool in sorted(_registered().items()):
        fn = _fn(tool)
        hints = typing.get_type_hints(fn)
        for parameter in inspect.signature(fn).parameters:
            hint = hints[parameter]
            options = (hint, *typing.get_args(hint))
            if any(typing.get_origin(option) is list for option in options):
                found.append((name, parameter))
    return found


@pytest.mark.parametrize("tool_name,parameter", collection_parameters())
@pytest.mark.parametrize("values", [
    ["x"] * (MAX_COLLECTION_ITEMS + 1),
    ["x" * (MAX_COLLECTION_ITEM_CHARS + 1)],
    ["x" * MAX_COLLECTION_ITEM_CHARS] * (MAX_COLLECTION_CHARS // MAX_COLLECTION_ITEM_CHARS + 1),
])
def test_every_collection_rejects_before_index_lookup(monkeypatch, tool_name, parameter, values):
    # Without this gate _need_index raises, or a corpus-dependent operation
    # would run. Exercise actual tool bodies, not only the shared helper.
    monkeypatch.setattr(app, "_INDEX", None)
    fn = _fn(_registered()[tool_name])
    kwargs = {name: "placeholder" for name, p in inspect.signature(fn).parameters.items()
              if p.default is inspect.Parameter.empty}
    for name in kwargs:
        if (tool_name, name) in collection_parameters():
            kwargs[name] = ["placeholder"]
    kwargs[parameter] = values
    result = fn(**kwargs)
    row = result[0] if isinstance(result, list) else result
    assert row["code"] == "invalid_argument"
    assert parameter in row["error"]
    assert set(row) == {"error", "code"}


def test_collection_inventory_is_complete():
    assert len(collection_parameters()) == 9


def test_exact_budgets_omission_and_empty_are_not_truncated():
    for values in (None, [], ["x"] * MAX_COLLECTION_ITEMS,
                   ["x" * MAX_COLLECTION_ITEM_CHARS] * (MAX_COLLECTION_CHARS // MAX_COLLECTION_ITEM_CHARS)):
        before = None if values is None else list(values)
        assert _validate_collection(values, "selector") is None
        assert values == before


@pytest.mark.parametrize("values", ["not-a-list", {"x": "y"}, [42]])
def test_invalid_collection_types(values):
    with pytest.raises(ValueError):
        _validate_collection(values, "selector")


def test_maximum_batch_preserves_order_and_duplicates(monkeypatch):
    from types import SimpleNamespace
    from mcpsrv.tools.papers import get_papers
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers={}))
    hashes = [str(n % 3) for n in range(MAX_COLLECTION_ITEMS)]
    assert [row["hash"] for row in get_papers(hashes)] == hashes


@pytest.mark.parametrize("fields", [[], ["unknown"], ["title", "first_author", "year"]])
def test_paper_header_projection_does_not_read_document_artifacts(monkeypatch, fields):
    from pathlib import Path
    from types import SimpleNamespace
    from mcpsrv.tools import papers
    record = {"hash": "abc", "hash_dir": "/unused", "title": "Title", "year": 2000,
              "authors": [{"surname": "Author"}]}
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers={"abc": record}))
    def forbidden(*args, **kwargs):
        pytest.fail("Header-only projection must not read or inventory document files")
    monkeypatch.setattr(papers, "_load_json", forbidden)
    monkeypatch.setattr(Path, "glob", forbidden)
    expected = {k: v for k, v in {**record, "first_author": "Author"}.items() if k in fields}
    expected["hash"] = "abc"
    assert papers.get_papers(["abc", "missing", "abc"], fields=fields) == [
        expected, {"hash": "missing", "error": "not_found"}, expected]


@pytest.mark.parametrize("field,expected_reads", [("top_taxa", ["taxa.json"]),
                                                 ("top_lexicon_terms", ["anatomy.json"])])
def test_paper_detail_projection_reads_only_its_requested_artifacts(tmp_path, monkeypatch, field, expected_reads):
    import json
    from types import SimpleNamespace
    from mcpsrv.tools import papers
    (tmp_path / "taxa.json").write_text(json.dumps({"taxa": [{"accepted_name": "Taxon"}]}))
    (tmp_path / "anatomy.json").write_text(json.dumps({"category": "anatomy", "terms": [{"canonical": "term"}]}))
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers={"abc": {
        "hash": "abc", "hash_dir": str(tmp_path), "title": "Title", "authors": []}}))
    full = papers.get_papers(["abc"])[0]
    original = papers._load_json
    reads = []
    def read(path, **kwargs):
        reads.append(path.name)
        return original(path, **kwargs)
    monkeypatch.setattr(papers, "_load_json", read)
    assert papers.get_papers(["abc"], fields=[field]) == [{"hash": "abc", field: full[field]}]
    assert reads == expected_reads

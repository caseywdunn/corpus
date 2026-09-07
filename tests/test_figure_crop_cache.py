"""All figure paths crop outside the bundle and share bounded, fresh results."""
import asyncio
from concurrent.futures import ThreadPoolExecutor
import io
import json

import pytest
from PIL import Image

from mcpsrv.figure_cache import FigureCropCache, figure_path
from mcpsrv.figure_http import make_figure_app
from mcpsrv.tools.figures import get_figure_image, get_figure_roi_image
from tests.test_figure_licensing_states import _make_index, _CLEARED, HASH


def fetch(app, query=b"profile=report&label=A"):
    messages = []
    async def receive():
        return {"type": "http.request", "body": b""}
    async def send(message):
        messages.append(message)
    asyncio.run(app({"type": "http", "method": "GET", "path": f"/figures/{HASH}/docling_1",
                     "query_string": query, "headers": []}, receive, send))
    return messages[0]["status"], b"".join(m.get("body", b"") for m in messages[1:])


def test_http_crops_on_first_request_with_readonly_bundle(tmp_path):
    idx = _make_index(tmp_path, work=_CLEARED)
    paths = list(tmp_path.rglob("*"))
    before = {p: p.read_bytes() for p in paths if p.is_file()}
    for p in paths:
        p.chmod(0o555 if p.is_dir() else 0o444)
    tmp_path.chmod(0o555)
    try:
        status, data = fetch(make_figure_app(idx))
        assert status == 200
        assert Image.open(io.BytesIO(data)).size == (20, 20)
        info = get_figure_roi_image(HASH, "docling_1", "A")
        assert info["crop"] is True
        image = get_figure_image(HASH, "docling_1", "A")
        assert image.data == data
        assert get_figure_image(HASH, "docling_1").path is not None
        assert {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()} == before
        assert not (tmp_path / "documents" / HASH / "figures" / "crops").exists()
    finally:
        tmp_path.chmod(0o755)
        for p in paths:
            p.chmod(0o755 if p.is_dir() else 0o644)


def test_cache_changes_with_pixels_and_roi_even_when_filename_same(tmp_path):
    src = tmp_path / "source.png"
    Image.new("RGB", (20, 20), "red").save(src)
    cache = FigureCropCache()
    first, a = cache.crop(src, [0, 0, 10, 10])
    assert cache.crop(src, [0, 0, 10, 10]) == (first, a)
    other_roi, b = cache.crop(src, [0, 0, 5, 5])
    assert other_roi != first and a != b
    Image.new("RGB", (20, 20), "blue").save(src)
    new_source, c = cache.crop(src, [0, 0, 10, 10])
    assert new_source != first and a != c


def test_concurrent_requests_and_eviction_are_bounded(tmp_path):
    src = tmp_path / "source.png"
    Image.new("RGB", (20, 20), "red").save(src)
    cache = FigureCropCache(max_entries=2, max_bytes=1000)
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda _: cache.crop(src, [0, 0, 10, 10]), range(12)))
    assert len({p for p, _ in results}) == 1
    first = results[0][0]
    for n in (5, 8, 12):
        cache.crop(src, [0, 0, n, n])
    assert not first.exists()
    assert len(list(cache.root.iterdir())) == 2
    assert sum(p.stat().st_size for p in cache.root.iterdir()) <= 1000
    assert cache.crop(src, [0, 0, 10, 10])[1] == results[0][1]


def test_unknown_roi_falls_back_to_whole_without_a_cache(tmp_path):
    idx = _make_index(tmp_path, work=_CLEARED)
    status, data = fetch(make_figure_app(idx), b"profile=report&label=unknown")
    assert status == 200 and Image.open(io.BytesIO(data)).size == (40, 40)
    assert not hasattr(idx, "figure_crop_cache")


def test_traversal_in_figure_record_cannot_read_another_file(tmp_path):
    idx = _make_index(tmp_path, work=_CLEARED)
    hd = tmp_path / "documents" / HASH
    payload = json.loads((hd / "figures.json").read_text())
    payload["figures"][0]["filename"] = "../figures.json"
    (hd / "figures.json").write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="inside"):
        figure_path(hd, payload["figures"][0])
    assert get_figure_roi_image(HASH, "docling_1", "A")["code"] == "unavailable"
    assert fetch(make_figure_app(idx))[0] == 400


def test_oversized_crop_is_rejected_before_cache_write(tmp_path, monkeypatch):
    from mcpsrv import figure_cache
    src = tmp_path / "source.png"
    Image.new("RGB", (20, 20), "red").save(src)
    cache = FigureCropCache()
    monkeypatch.setattr(figure_cache, "MAX_IMAGE_PIXELS", 100)
    with pytest.raises(ValueError, match="pixel limit"):
        cache.crop(src, [0, 0, 10, 10])
    assert not list(cache.root.iterdir())

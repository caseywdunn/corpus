"""Bounded disposable figure crops; never write into an immutable bundle."""
from collections import OrderedDict
import hashlib
import io
import json
from pathlib import Path
import tempfile
import threading

MAX_SOURCE_BYTES = 128 * 1024 * 1024
MAX_IMAGE_PIXELS = 64_000_000
_INITIALIZE_LOCK = threading.Lock()


def figure_path(hash_dir, figure):
    root = (Path(hash_dir) / "figures").resolve(strict=True)
    path = (root / (figure.get("filename") or "")).resolve(strict=True)
    if not path.is_file() or not path.is_relative_to(root):
        raise ValueError("figure path is not a file inside its figure directory")
    return path


class FigureCropCache:
    """Process-owned LRU with a byte/count ceiling and content+ROI keys."""

    def __init__(self, *, max_bytes=128 * 1024 * 1024, max_entries=128):
        self._temporary = tempfile.TemporaryDirectory(prefix="corpus-figure-crops-")
        self.root = Path(self._temporary.name)
        self.max_bytes = max_bytes
        self.max_entries = max_entries
        self._entries = OrderedDict()
        self._bytes = 0
        self._lock = threading.Lock()

    def crop(self, source, roi):
        from PIL import Image
        if not isinstance(roi, (list, tuple)) or len(roi) != 4:
            raise ValueError("figure ROI must contain four coordinates")
        coordinates = tuple(int(v) for v in roi)
        with self._lock:
            with Path(source).open("rb") as stream:
                raw = stream.read(MAX_SOURCE_BYTES + 1)
            if len(raw) > MAX_SOURCE_BYTES:
                raise ValueError("figure exceeds the server source-image byte limit")
            key = hashlib.sha256(raw + json.dumps(coordinates).encode()).hexdigest()
            path = self.root / f"crop-{key}.png"
            if key in self._entries and path.is_file():
                self._entries.move_to_end(key)
                return path, path.read_bytes()
            if key in self._entries:
                self._bytes -= self._entries.pop(key)
            with Image.open(io.BytesIO(raw)) as image:
                if image.width * image.height > MAX_IMAGE_PIXELS:
                    raise ValueError("figure exceeds the server source-image pixel limit")
                x0, y0, x1, y1 = coordinates
                x0 = max(0, min(x0, image.width - 1))
                y0 = max(0, min(y0, image.height - 1))
                x1 = max(x0 + 1, min(x1, image.width))
                y1 = max(y0 + 1, min(y1, image.height))
                buffer = io.BytesIO()
                image.crop((x0, y0, x1, y1)).save(buffer, format="PNG")
            data = buffer.getvalue()
            if len(data) > self.max_bytes or self.max_entries < 1:
                raise ValueError("crop exceeds disposable cache limits")
            while self._entries and (len(self._entries) >= self.max_entries or self._bytes + len(data) > self.max_bytes):
                stale, size = self._entries.popitem(last=False)
                (self.root / f"crop-{stale}.png").unlink(missing_ok=True)
                self._bytes -= size
            path.write_bytes(data)  # Private tree; every reader holds this lock.
            self._entries[key] = len(data)
            self._bytes += len(data)
            return path, data


def crop_figure(idx, source, roi):
    with _INITIALIZE_LOCK:
        cache = getattr(idx, "figure_crop_cache", None)
        if cache is None:
            cache = FigureCropCache()
            idx.figure_crop_cache = cache
    return cache.crop(source, roi)

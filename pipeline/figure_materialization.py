"""Reset derived figures from the prepared PDF, never from prior split state.

Panel-mode/model changes are rare build operations. Re-extract into a scratch
directory before replacing the figure base, avoiding a permanent second image
cache. Figure-only resets do not rewrite stored text/chunks or rerun OCR and
bibliographic extraction. Full extraction uses the same publication path.
"""

import json
from pathlib import Path
import tempfile


def rebuild_figure_base(hash_dir: Path, extract, *, figures_only=True):
    """Install freshly extracted figures, rolling back on ordinary I/O failure.

    With ``figures_only=False``, publish text and the optional Docling sidecar
    too. This is the full extraction path: no stale images or sidecar from a
    previous generation survive. Call inside a cleared producer receipt.
    If killed during publication, retry reconstructs the base from processed.pdf; it must
    not treat whatever split/renamed files remain as source evidence.
    """
    pdf = hash_dir / "processed.pdf"
    if not pdf.is_file():
        raise FileNotFoundError(f"{pdf}: re-extract the document before resetting figures")
    scan_path = hash_dir / "scan_detection.json"
    scan = json.loads(scan_path.read_text()) if scan_path.exists() else {}
    figures_dir = hash_dir / "figures"
    with tempfile.TemporaryDirectory(prefix=".figure-reset-", dir=hash_dir) as tmp:
        scratch = Path(tmp)
        fresh_dir = scratch / "figures"
        fresh_dir.mkdir()
        fresh_json = scratch / "figures.json"
        extract(pdf, scratch / "text.json", fresh_json, fresh_dir,
                scan_file_type=scan.get("file_type"))
        data = json.loads(fresh_json.read_text())
        if not isinstance(data.get("figures"), list):
            raise ValueError("Fresh extraction did not produce a figures list")
        for figure in data["figures"]:
            filename = figure.get("filename")
            if not filename or Path(filename).name != filename or filename in (".", ".."):
                raise ValueError(f"Unsafe extracted figure filename: {filename!r}")
            if not (fresh_dir / filename).is_file():
                raise ValueError(f"Extracted figure image is missing: {filename}")
            figure["file_path"] = str(figures_dir / filename)
        data["figures_directory"] = str(figures_dir)
        fresh_json.write_text(json.dumps(data, indent=2))
        names = ["figures", "figures.json"]
        if not figures_only:
            if not (scratch / "text.json").is_file():
                raise ValueError("Fresh extraction did not produce text.json")
            names += ["text.json", "docling_doc.json"]
        backups = []
        installed = []
        try:
            for name in names:
                destination = hash_dir / name
                if destination.exists():
                    destination.rename(scratch / f"previous-{name}")
                    backups.append(name)
                fresh = scratch / name
                if fresh.exists():
                    fresh.rename(destination)
                    installed.append(name)
        except BaseException:
            for name in reversed(installed):
                (hash_dir / name).rename(scratch / f"failed-{name}")
            for name in reversed(backups):
                (scratch / f"previous-{name}").rename(hash_dir / name)
            raise

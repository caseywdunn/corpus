"""Bounded source-gap corroboration for missing spaces in native PDF text.

Geometry proposes boundaries, two independent segmentation modes corroborate
only those boundaries. OCR never supplies letters or boundaries on its own.
This build-time fallback leaves uncertain candidates unchanged and reviewable.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import importlib.metadata
import os
from pathlib import Path
import re
import shutil
from statistics import median
import subprocess

SPACE_POLICY = "exact-letters-geometric-gaps-ocr-intersection-v2"


@lru_cache(maxsize=4)
def _ocr_installation(executable, modified, tessdata_prefix):
    del modified, tessdata_prefix
    try:
        version = subprocess.run([executable, "--version"], capture_output=True,
                                 text=True, timeout=5, check=True).stdout.strip()
        languages = subprocess.run([executable, "--list-langs"], capture_output=True,
                                   text=True, timeout=5, check=True).stdout
        location = re.search(r'"([^"]+)"', languages)
        model = Path(location[1]) / "eng.traineddata" if location else None
        if model is None or not model.is_file():
            return version, None
        return version, str(model)
    except (OSError, subprocess.SubprocessError):
        return None, None


@lru_cache(maxsize=4)
def _model_hash(path, modified, size):
    del modified, size
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def source_spacing_producer():
    """Cheap offline producer receipt; no OCR, model loading, or network."""
    result = {"policy": SPACE_POLICY, "ocr_modes": [6, 7], "dpi": 400,
              "language": "eng", "max_candidates_per_page": 8,
              "timeout_seconds": 15, "available": False,
              "pymupdf_version": importlib.metadata.version("pymupdf")}
    executable = shutil.which("tesseract")
    if not executable:
        return result
    try:
        version, model = _ocr_installation(executable, Path(executable).stat().st_mtime_ns,
                                           os.environ.get("TESSDATA_PREFIX"))
        if model is None:
            return {**result, "version": version}
        stat = Path(model).stat()
        return {**result, "available": True, "version": version,
                "traineddata_sha256": _model_hash(model, stat.st_mtime_ns, stat.st_size)}
    except OSError:
        return result


def geometric_space_candidate(line):
    """Return native lines with strongly separated, uniform gap clusters.

    Ordinary letter spacing, kerning, font changes, and uniform OCR character
    boxes provide no such evidence. The thresholds are deliberately conservative;
    a missed repair is retained for review rather than introducing a word split.
    """
    if tuple(line.get("dir", (1, 0))) != (1, 0):
        return None
    chars = []
    for span in line.get("spans", []):
        chars.extend((char, span.get("font"), span.get("size")) for char in span.get("chars", []))
    text = "".join(c[0]["c"] for c in chars)
    choices = []
    isolated = False
    for match in re.finditer(r"[^\W\d_]{20,}", text):
        # A PDF can omit the space between an italic name and upright prose.
        # Evaluate existing uniform-font portions of that joined native run;
        # a font transition itself is never evidence for an inserted space.
        start = match.start()
        portions = []
        for end in range(start + 1, match.end()):
            if chars[end][1:] != chars[end - 1][1:]:
                portions.append((start, end))
                start = end
        portions.append((start, match.end()))
        for start, end in portions:
            if end - start < 20:
                continue
            choice = _geometric_choice(chars, text, start, end)
            if choice:
                choices.append(choice)
                isolated |= len(portions) > 1
    if choices:
        candidate = {"text": text, "bbox": list(line["bbox"]), "choices": choices}
        if isolated:
            # One crop per source line still covers every proposed run. Avoid
            # making an unrelated style transition part of OCR segmentation.
            boxes = [char["bbox"] for choice in choices
                     for char, _, _ in chars[choice["start"]:choice["end"]]]
            candidate["ocr_bbox"] = [min(b[0] for b in boxes), min(b[1] for b in boxes),
                                     max(b[2] for b in boxes), max(b[3] for b in boxes)]
            candidate["isolated_font_runs"] = True
        return candidate
    return None


def _geometric_choice(chars, text, start, end):
    """Measure only internal gaps of one existing uniform-font letter run."""
    subset = chars[start:end]
    widths = [c["bbox"][2] - c["bbox"][0] for c, _, _ in subset]
    width = median(widths)
    if width <= 0:
        return None
    gaps = [right[0]["bbox"][0] - left[0]["bbox"][2]
            for left, right in zip(subset, subset[1:])]
    baseline = median(gaps)
    if baseline > width * .06 or baseline < -width * .06:
        return None
    threshold = max(.65, width * .2, baseline * 4 + .3)
    boundaries = [i + 1 for i, gap in enumerate(gaps) if gap >= threshold]
    high = [gaps[i - 1] for i in boundaries]
    if not 3 <= len(high) <= len(gaps) * .3 or max(high) > min(high) * 1.25:
        return None
    candidate = text[start:end]
    for i in reversed(boundaries):
        candidate = candidate[:i] + " " + candidate[i:]
    return {"original": text[start:end], "candidate": candidate,
            "start": start, "end": end,
            "boundaries": boundaries, "gaps": [round(v, 4) for v in high],
            "median_character_width": round(width, 4),
            "median_intra_word_gap": round(baseline, 4)}


def _supported_boundaries(original, ocr_text):
    compact = "".join(c for c in ocr_text if not c.isspace())
    positions = [i for i, c in enumerate(ocr_text) if not c.isspace()]
    start = compact.find(original)
    if start < 0 or compact.find(original, start + 1) >= 0:
        return set()
    return {i for i in range(1, len(original))
            if any(c.isspace() for c in ocr_text[positions[start+i-1]+1:positions[start+i]])}


def corroborate_source_gaps(page, candidates, producer):
    """Return verified extra source lines and bounded failure observations."""
    import fitz

    lines = []
    observations = []
    for index, candidate in enumerate(candidates):
        observation = {**candidate, "route": SPACE_POLICY}
        if index >= producer["max_candidates_per_page"]:
            observations.append({**observation, "status": "candidate_budget_exceeded"})
            continue
        if not producer["available"]:
            observations.append({**observation, "status": "ocr_unavailable"})
            continue
        rect = fitz.Rect(candidate.get("ocr_bbox", candidate["bbox"])) + (-2, -2, 2, 2)
        rect &= page.rect
        # No page raster retained; only one small candidate line at a time.
        if rect.width * rect.height * (producer["dpi"] / 72) ** 2 > 4_000_000:
            observations.append({**observation, "status": "crop_pixel_budget_exceeded"})
            continue
        png = page.get_pixmap(clip=rect, dpi=producer["dpi"], alpha=False).tobytes("png")
        observation.update({"crop_bbox": list(rect), "crop_dpi": producer["dpi"],
                            "crop_sha256": hashlib.sha256(png).hexdigest()})
        outputs = []
        try:
            for mode in producer["ocr_modes"]:
                result = subprocess.run([shutil.which("tesseract"), "stdin", "stdout",
                                         "--psm", str(mode), "-l", producer["language"]],
                                        input=png, capture_output=True, check=True,
                                        timeout=producer["timeout_seconds"])
                outputs.append(result.stdout.decode("utf-8").strip())
        except (OSError, subprocess.SubprocessError) as exc:
            observations.append({**observation, "status": "ocr_failed", "error": type(exc).__name__})
            continue
        repaired = candidate["text"]
        accepted = []
        for choice in reversed(candidate["choices"]):
            wanted = set(choice["boundaries"])
            # A complete candidate's geometric boundaries must all have both
            # OCR confirmations. OCR-only splits (e.g. 'cou nt') are ignored.
            if wanted and all(wanted <= _supported_boundaries(choice["original"], output)
                              for output in outputs):
                repaired = repaired[:choice["start"]] + choice["candidate"] + repaired[choice["end"]:]
                accepted.append(choice)
        observations.append({**observation, "ocr_outputs": outputs,
                             "accepted": accepted, "status": "verified" if accepted else "ocr_disagreement"})
        if accepted:
            lines.append(repaired)
    return lines, observations

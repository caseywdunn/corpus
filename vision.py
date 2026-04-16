"""Vision-model backends for Pass 3b figure panel / sub-figure detection.

PLAN.md §9. Pass 3a (Tesseract OCR in figures.py) has low recall on
line-art scientific figures — around 20–40% of panels on the demo set.
Pass 3b escalates to a vision-language model that can read embedded
labels (``A``, ``B``, ``C``) and detect compound figures (``Fig. 3`` +
``Fig. 4`` in one image) much more reliably.

Backend abstraction mirrors ``embeddings.py``:

* :class:`ClaudeVisionBackend` — Anthropic Claude via the official SDK.
  Network-dependent, pennies per figure at Haiku. The right choice for
  development and any host with outbound HTTP to the Anthropic API.
* :class:`LocalVLMBackend` (forthcoming) — Qwen2.5-VL or similar
  open-weights model on CUDA / MPS / CPU. Zero per-call cost, fully
  local, network-independent. The right choice for the Bouchet
  production run.

Both backends return the same structured output — a list of panel /
embedded-figure ROIs with normalized bboxes + per-ROI confidence — so
the Pass 3b pipeline step is backend-agnostic.
"""

from __future__ import annotations

import base64
import io
import json
import logging
import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional

from dotenv import load_dotenv
load_dotenv()  # picks up ANTHROPIC_API_KEY from .env at import time

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Output contract (what every backend must return)
# ---------------------------------------------------------------------------
#
# detect_figure_panels() returns a list of dicts:
#
#   {
#     "type": "panel" | "embedded_figure",
#     "label": "A" | ... (for panels),
#     "figure_number": "4" | ... (for embedded_figure),
#     "bbox_px": [x0, y0, x1, y1],    # top-left origin, PIL convention
#     "confidence": float,             # 0..1
#     "description": str,              # brief description of this panel
#     "source": "vision:claude-haiku-4-5" | "vision:qwen2.5-vl-7b",
#   }
#
# Pass 3b in figures.py reconciles these with panels_from_caption so the
# final figures.json entry looks the same whether Pass 3a or Pass 3b
# populated it.


class VisionBackendError(RuntimeError):
    """Raised when the vision backend fails irrecoverably (bad API key,
    model not loadable, non-retryable API error). The caller typically
    records ``pass3_status = "vision_backend_failed"`` and moves on."""


class VisionBackend(ABC):
    """Common interface every Pass 3b backend implements."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Short identifier recorded on each emitted ROI's ``source`` field
        (e.g., ``"vision:claude-haiku-4-5"``, ``"vision:qwen2.5-vl-7b"``).
        Used downstream to tell where an annotation came from."""

    @abstractmethod
    def detect_figure_panels(
        self,
        image_path: Path,
        caption_text: str,
        expected_labels: List[str],
    ) -> List[Dict]:
        """Return the panel / embedded-figure ROI list for one image."""


# ---------------------------------------------------------------------------
# Claude backend
# ---------------------------------------------------------------------------


# Kept as a module-level constant so the prompt is easy to review without
# spelunking through methods. The prompt is terse on purpose — Claude is
# already excellent at this task and responds better to clear constraints
# than to long elaborations.
_CLAUDE_SYSTEM_PROMPT = """You are analyzing a scientific figure extracted from a biology paper.

Your task has two parts:

  (1) Panel labels — single capital letters (A, B, C, ...) that label sub-panels
      within a single figure.
  (2) Compound detection — when the image actually contains two or more separate
      figures merged into one extracted image (e.g., Fig. 3 and Fig. 4 side by
      side). Signals: the number of panel labels is LARGER than the caption's
      expected panel count, or duplicate panel letters appear in the image (two
      separate As, two Bs, etc.), or visible "Fig. N" text labels with different
      N values.

Return STRICTLY VALID JSON:

{
  "panels": [
    {
      "label": "A",
      "parent_figure_index": 0,         // 0 = the first sub-figure, 1 = the second, etc.
                                         // Use 0 for every panel when the image is not compound.
      "panel_bbox_norm": [x0, y0, x1, y1],   // the WHOLE panel's bounding box
      "label_bbox_norm": [x0, y0, x1, y1],   // the small bbox of just the letter label
      "confidence": 0.0-1.0,             // your confidence the label is visible and placed correctly
      "description": "short phrase describing the panel"
    }
  ],
  "embedded_figures": [
    {
      "parent_figure_index": 1,         // matches the parent_figure_index on panels above
      "figure_number": "4",              // the N from a "Fig. N" label if visible; null if inferred
      "panel_bbox_norm": [x0, y0, x1, y1],   // bbox of the entire sub-figure region
      "confidence": 0.0-1.0
    }
  ]
}

Coordinate rules:
- All bboxes are [left, top, right, bottom] in NORMALIZED image coordinates
  (each coordinate is a float in 0.0 .. 1.0). Top-left origin (PIL convention,
  y grows downward).
- panel_bbox_norm must cover the ENTIRE panel that the label identifies — the
  image region that would remain if you cropped the panel out of the figure.
  Typically panel labels sit at a corner of their panel; your panel_bbox_norm
  should extend from that corner outward to cover the panel's content.
- label_bbox_norm is the small bbox of just the letter itself.

Interpretation rules:
- ONLY emit a panel when the letter label is actually visible in the image.
  Do not infer panels from the caption alone.
- If the caption expects panels A-B (2) but you see labels A, B, A, B in the
  image, this is a compound of two separate figures. Assign
  parent_figure_index = 0 to the first {A, B} set and parent_figure_index = 1
  to the second; emit one embedded_figures entry per sub-figure.
- If the image has no recognizable panel or figure labels, return
  {"panels": [], "embedded_figures": []}.
- Output ONLY the JSON. No preface, no markdown fence, no explanation."""


def _extract_json(text: str) -> Optional[dict]:
    """Pull the outermost JSON object from an LLM response.

    Claude usually obeys "output ONLY JSON" but sometimes wraps it in
    markdown fences or adds prose. Strip ``` fences and take the first
    balanced ``{ ... }`` span.
    """
    if not text:
        return None
    # Strip markdown fences.
    text = re.sub(r"^```(?:json)?\s*", "", text.strip(), flags=re.MULTILINE)
    text = re.sub(r"\s*```$", "", text.strip(), flags=re.MULTILINE)
    # Find first balanced {...}.
    depth = 0
    start = None
    for i, c in enumerate(text):
        if c == "{":
            if depth == 0:
                start = i
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0 and start is not None:
                try:
                    return json.loads(text[start:i + 1])
                except json.JSONDecodeError:
                    return None
    return None


class ClaudeVisionBackend(VisionBackend):
    """Pass 3b backend using Anthropic Claude's vision-enabled Messages API.

    ``model`` defaults to Haiku 4.5 — fastest and cheapest of the Claude
    4-family models, accurate enough for the panel-detection task in
    all the demo spot-checks. Bump to a Sonnet / Opus model if Haiku's
    quality falls short on some subset of figures.

    Requires ``ANTHROPIC_API_KEY`` in the environment (or in ``.env``,
    which is loaded at module import).
    """

    def __init__(
        self,
        model: str = "claude-haiku-4-5-20251001",
        max_tokens: int = 1024,
    ):
        if not os.environ.get("ANTHROPIC_API_KEY"):
            raise VisionBackendError(
                "ANTHROPIC_API_KEY not set. Put it in .env at the repo root or "
                "export it in your shell."
            )
        try:
            import anthropic
        except ImportError as e:
            raise VisionBackendError(
                "anthropic package not installed (pip install anthropic)"
            ) from e
        self._anthropic = anthropic
        self.client = anthropic.Anthropic()
        self.model = model
        self.max_tokens = max_tokens

    @property
    def name(self) -> str:
        # Strip the date suffix from the model id for cleaner ROI sources.
        stem = self.model.rsplit("-", 1)[0] if self.model.count("-") >= 3 else self.model
        return f"vision:{stem}"

    def _encode_image(self, image_path: Path) -> tuple:
        """Load image, return (base64_png, width, height). Resizes if
        larger than ~2000 px max-side — Claude handles up to 8000 but
        paying for a high-resolution render of a small figure is wasteful
        and doesn't improve label detection."""
        from PIL import Image
        with Image.open(image_path) as img:
            img = img.convert("RGB")
            w, h = img.size
            max_side = max(w, h)
            if max_side > 2000:
                ratio = 2000 / max_side
                img = img.resize((int(w * ratio), int(h * ratio)))
                w, h = img.size
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            b64 = base64.standard_b64encode(buf.getvalue()).decode()
        return b64, w, h

    def detect_figure_panels(
        self,
        image_path: Path,
        caption_text: str,
        expected_labels: List[str],
    ) -> List[Dict]:
        try:
            image_b64, w, h = self._encode_image(image_path)
        except Exception as e:
            raise VisionBackendError(f"could not read image {image_path}: {e}") from e

        # User message includes the caption + expected labels so Claude
        # can ground its bbox hunt in what's supposed to be there.
        user_text = (
            f"Caption of this figure: {caption_text!r}\n\n"
            f"Expected panel labels from the caption (confirm visibility before "
            f"emitting): {expected_labels}\n\n"
            f"Image dimensions (px): {w} × {h}."
        )

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                system=_CLAUDE_SYSTEM_PROMPT,
                messages=[{
                    "role": "user",
                    "content": [
                        {
                            "type": "image",
                            "source": {
                                "type": "base64",
                                "media_type": "image/png",
                                "data": image_b64,
                            },
                        },
                        {"type": "text", "text": user_text},
                    ],
                }],
            )
        except self._anthropic.APIError as e:
            raise VisionBackendError(f"Claude API error: {e}") from e

        # Concatenate all text blocks; usually there's just one.
        response_text = "".join(
            b.text for b in response.content if getattr(b, "type", None) == "text"
        )
        parsed = _extract_json(response_text)
        if not parsed:
            logger.warning(
                "Could not extract JSON from Claude response for %s; "
                "response head: %r",
                image_path.name, response_text[:200],
            )
            return []

        def _norm_to_px(bbox_norm):
            if not (isinstance(bbox_norm, list) and len(bbox_norm) == 4):
                return None
            try:
                x0 = int(max(0.0, float(bbox_norm[0])) * w)
                y0 = int(max(0.0, float(bbox_norm[1])) * h)
                x1 = int(min(1.0, float(bbox_norm[2])) * w)
                y1 = int(min(1.0, float(bbox_norm[3])) * h)
            except (TypeError, ValueError):
                return None
            if x1 <= x0 or y1 <= y0:
                return None
            return [x0, y0, x1, y1]

        out: List[Dict] = []
        src = self.name
        for p in parsed.get("panels") or []:
            panel_px = _norm_to_px(p.get("panel_bbox_norm"))
            if panel_px is None:
                continue
            entry = {
                "type": "panel",
                "label": p.get("label", ""),
                "parent_figure_index": int(p.get("parent_figure_index", 0)),
                "bbox_px": panel_px,
                "confidence": float(p.get("confidence", 0.0)),
                "description": (p.get("description") or "").strip(),
                "source": src,
            }
            label_px = _norm_to_px(p.get("label_bbox_norm"))
            if label_px is not None:
                entry["label_bbox_px"] = label_px
            out.append(entry)
        for f in parsed.get("embedded_figures") or []:
            panel_px = _norm_to_px(f.get("panel_bbox_norm"))
            if panel_px is None:
                continue
            out.append({
                "type": "embedded_figure",
                "parent_figure_index": int(f.get("parent_figure_index", 0)),
                "figure_number": str(f.get("figure_number") or "").strip() or None,
                "bbox_px": panel_px,
                "confidence": float(f.get("confidence", 0.0)),
                "source": src,
            })
        return out


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def get_vision_backend(backend: str = "claude", **kwargs) -> VisionBackend:
    """Construct a vision backend by name. ``"claude"`` for the Anthropic
    API backend. ``"local"`` reserved for a future open-weights backend
    (Qwen2.5-VL on CUDA/MPS)."""
    backend = (backend or "claude").lower()
    if backend == "claude":
        return ClaudeVisionBackend(**kwargs)
    if backend == "local":
        raise VisionBackendError(
            "LocalVLMBackend not yet implemented; use --vision-backend claude for now"
        )
    raise ValueError(f"Unknown vision backend: {backend!r}")

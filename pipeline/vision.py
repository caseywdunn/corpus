"""Vision-model backends for Pass 3b figure panel / sub-figure detection.

See dev_docs/OVERVIEW.md "Figure pipeline". Pass 3a (Tesseract OCR in
pipeline/figures.py) has low recall on
line-art scientific figures — around 20–40% of panels on the demo set.
Pass 3b escalates to a vision-language model that can read embedded
labels (``A``, ``B``, ``C``) and detect compound figures (``Fig. 3`` +
``Fig. 4`` in one image) much more reliably.

Backend abstraction mirrors ``pipeline.embeddings``:

* :class:`ClaudeVisionBackend` — Anthropic Claude via the official SDK.
  Network-dependent, pennies per figure at Haiku. The right choice for
  development and any host with outbound HTTP to the Anthropic API.
* :class:`LocalVLMBackend` — Qwen2.5-VL or similar open-weights model
  on CUDA / MPS / CPU. Zero per-call cost, fully local, network-
  independent. The right choice for the Bouchet production run.

Both backends return the same structured output — a list of panel /
embedded-figure ROIs with normalized bboxes + per-ROI confidence — so
the Pass 3b pipeline step is backend-agnostic.
"""

from __future__ import annotations

import base64
import io
import json
from collections import Counter
from functools import cached_property
import logging
import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional

from dotenv import load_dotenv
from .model_provenance import DEFAULT_VISION_MODELS, vision_producer
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
# Pass 3b in pipeline/figures.py reconciles these with panels_from_caption so the
# final figures.json entry looks the same whether Pass 3a or Pass 3b
# populated it.


class VisionBackendError(RuntimeError):
    """Raised when a backend cannot return a complete usable answer.

    This includes setup/API failures, malformed output, and token-budget
    truncation. The caller records ``pass3_status = "vision_backend_failed"``
    rather than confusing any of those with a clean no-label result.
    """


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
        """Return the panel / embedded-figure ROI list for one image.

        ``expected_labels`` is normally a list of panel letters. For a known
        grouped plate it is instead a list of caption-enumerated figure
        numbers; those belong in ``embedded_figures``, not ``panels``. An
        empty list is reserved for explicitly admitted bare-plate discovery;
        the backend must emit only numbers it can actually see.
        """


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
  (2) Figure regions — when the image contains two or more separately numbered
      figures (e.g., Fig. 3 and Fig. 4 side by side), return each region under
      embedded_figures. This includes known grouped plates whose expected list
      contains figure numbers instead of panel letters. Other signals are more
      panel labels than the caption expects, duplicate panel letters (two As,
      two Bs), or visible figure-number labels.

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
- When the expected list contains numbers, look for those visible numbers and
  emit their whole numbered regions under embedded_figures. Do not reinterpret
  the numbers as panel letters and do not infer an unseen region.
- When the expected list is empty, the caller has explicitly admitted a bare
  historical plate. Discover Arabic numbers that visibly label distinct
  engravings and emit their whole regions under embedded_figures. Do not emit
  the plate number, page number, scale values, anatomical-key numbers,
  handwritten marks, or a number that is not visibly attached to a distinct
  engraving.
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
    balanced ``{ ... }`` span that actually parses.

    "That actually parses" is the part #253 needed: this used to return
    ``None`` on the first :class:`json.JSONDecodeError` and stop, so one
    malformed leading object discarded every well-formed one after it. A
    truncated response has no balanced span at all and still yields
    ``None`` — completion metadata is checked separately by the caller.
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
                    start = None          # keep scanning for a later span
    return None


# Output-budget model for the panel prompt (#253 / #269). One panel object
# carries six fields and measured 60–90 tokens; 120 is that worst case plus
# margin. The base covers the enclosing object and any unwanted preamble.
# Both backends use the same calculation: configured values are floors, not
# ceilings independent of the requested structure.
_VLM_TOKENS_PER_PANEL = 120
_VLM_TOKEN_BASE = 256
_VLM_DISCOVERY_TOKEN_FLOOR = 4096


def _vision_token_budget(floor: int, expected_labels) -> int:
    if not expected_labels:
        # Discovery does not know the output cardinality in advance. Historical
        # plates in the measured set carry up to the mid-teens of engravings;
        # a larger floor prevents valid JSON after an arbitrary prefix from
        # masquerading as a complete inventory.
        return max(int(floor), _VLM_DISCOVERY_TOKEN_FLOOR)
    return max(int(floor), _VLM_TOKEN_BASE + _VLM_TOKENS_PER_PANEL * len(
        expected_labels or []
    ))


def _vision_user_text(
    caption_text: str,
    expected_labels: List[str],
    width: int,
    height: int,
) -> str:
    """Compose the shared Claude/Qwen task-specific prompt."""
    if expected_labels:
        target_instruction = (
            "Expected panel letters OR grouped-plate figure numbers from the "
            "caption (confirm visibility before emitting): "
            f"{expected_labels}"
        )
    else:
        target_instruction = (
            "No figure-number list was recoverable from this bare plate "
            "caption. Discover every Arabic number visibly attached to a "
            "distinct engraving; emit the engraving's whole region under "
            "embedded_figures. Do not infer missing or unreadable numbers."
        )
    return (
        f"Caption of this figure: {caption_text!r}\n\n"
        f"{target_instruction}\n\n"
        f"Image dimensions (px): {width} × {height}."
    )


def _parse_complete_vision_response(
    response_text: str,
    *,
    backend: str,
    truncated: bool = False,
    truncation_detail: str = "",
) -> Dict:
    """Return a complete response payload or fail visibly (#269).

    A real empty detection is the explicit two-list object. Invalid JSON,
    a missing output collection, or a provider/model stop at its token cap is
    not evidence that the image contains no labels.
    """
    if truncated:
        detail = f" ({truncation_detail})" if truncation_detail else ""
        raise VisionBackendError(f"{backend} response was token-truncated{detail}")
    parsed = _extract_json(response_text)
    if not isinstance(parsed, dict):
        raise VisionBackendError(
            f"{backend} returned no complete JSON object; response head: "
            f"{response_text[:200]!r}"
        )
    for field in ("panels", "embedded_figures"):
        if not isinstance(parsed.get(field), list):
            raise VisionBackendError(
                f"{backend} returned an incomplete schema: {field!r} must be a list"
            )
    return parsed


# ---------------------------------------------------------------------------
# Model bboxes → pixels (#253)
# ---------------------------------------------------------------------------

# What units did the model actually emit? The prompt demands "each coordinate
# is a float in 0.0 .. 1.0", and both backends multiplied by the image
# dimensions on that assumption. Qwen2.5-VL frequently ignores it and emits
# absolute pixels instead — measured at 130 of 142 observable responses on
# one cluster, and 100% of those since 2026-05-30.
#
# The old conversion turned that into silent loss, two different ways:
#
#   [17, 808, 150, 1365] -> x0 = 17*w = 8432, x1 = min(1.0,150)*w = w
#                           x1 <= x0, so the panel was dropped, no warning.
#   [0, 0, 487, 606]     -> x0 = 0, x1 = w, y1 = h
#                           the guard passes and you get a "panel" spanning
#                           the whole figure — wrong data, not missing data,
#                           and it lands in the `completed` counter.
#
# A coordinate above 1.0 cannot be normalized, so the units are recoverable
# without guessing: use the values as pixels. Repairing strictly dominates
# dropping, and it removes the whole-figure impostor above.
_BBOX_NORMALIZED = "normalized"
_BBOX_PIXELS = "pixels"
_BBOX_OUT_OF_RANGE = "out_of_range"
_BBOX_MALFORMED = "malformed"
_BBOX_DEGENERATE = "degenerate"


def _bbox_to_px(bbox, w: int, h: int):
    """``(bbox_px | None, disposition)`` for one model-emitted bbox.

    ``disposition`` is one of the ``_BBOX_*`` constants and is counted by the
    caller, so a build reports what its model actually emitted instead of
    discarding the evidence. A bbox that survives is always ``[x0, y0, x1,
    y1]`` in pixels, clamped to the image.
    """
    if not (isinstance(bbox, list) and len(bbox) == 4):
        return None, _BBOX_MALFORMED
    try:
        vals = [float(v) for v in bbox]
    except (TypeError, ValueError):
        return None, _BBOX_MALFORMED

    if all(0.0 <= v <= 1.0 for v in vals):
        units = _BBOX_NORMALIZED
        x0, y0, x1, y1 = (vals[0] * w, vals[1] * h, vals[2] * w, vals[3] * h)
    elif all(0.0 <= v <= max(w, h) for v in vals):
        # Above 1.0 but inside the image: pixels. Note a bbox mixing the two
        # conventions reads as pixels and its sub-1.0 members collapse to 0 —
        # accepted, because a mixed bbox is not recoverable either way and
        # this at least keeps the panel.
        units = _BBOX_PIXELS
        x0, y0, x1, y1 = vals
    else:
        return None, _BBOX_OUT_OF_RANGE

    x0 = int(max(0, min(w, x0)))
    y0 = int(max(0, min(h, y0)))
    x1 = int(max(0, min(w, x1)))
    y1 = int(max(0, min(h, y1)))
    if x1 <= x0 or y1 <= y0:
        return None, _BBOX_DEGENERATE
    return [x0, y0, x1, y1], units


def _log_bbox_dispositions(counts, name: str, backend: str) -> None:
    """Say what the model emitted and what it cost, once per figure.

    The drop paths had no logging at all, which is why the units question
    could only be answered from the handful of responses that failed to
    parse — see #253. Pixel repairs log at info (they worked); anything lost
    logs at warning (it did not).
    """
    if not counts:
        return
    lost = sum(counts.get(k, 0) for k in
               (_BBOX_OUT_OF_RANGE, _BBOX_MALFORMED, _BBOX_DEGENERATE))
    repaired = counts.get(_BBOX_PIXELS, 0)
    if repaired:
        logger.info(
            "[%s] %s: %d bbox(es) arrived as pixels rather than the 0-1 "
            "floats the prompt asks for; converted. (%d normalized.)",
            name, backend, repaired, counts.get(_BBOX_NORMALIZED, 0),
        )
    if lost:
        logger.warning(
            "[%s] %s: dropped %d bbox(es) — %s. These become missing panel "
            "ROIs and the figure may be recorded as no-labels.",
            name, backend, lost,
            ", ".join(f"{k}={counts[k]}" for k in sorted(counts) if counts[k]),
        )


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
        model: str = DEFAULT_VISION_MODELS["vision-claude"],
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

    panel_mode = "vision-claude"

    @cached_property
    def producer(self):
        result = vision_producer(self.panel_mode, self.model)
        result["generation"]["max_tokens"] = self.max_tokens
        return result

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

    def _token_budget(self, expected_labels) -> int:
        """Scale the configured floor by the response structure requested."""
        return _vision_token_budget(self.max_tokens, expected_labels)

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
        user_text = _vision_user_text(caption_text, expected_labels, w, h)

        initial_budget = self._token_budget(expected_labels)
        for attempt, budget in enumerate(
            (initial_budget, initial_budget * 2), start=1,
        ):
            try:
                response = self.client.messages.create(
                    model=self.model,
                    max_tokens=budget,
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
                b.text for b in response.content
                if getattr(b, "type", None) == "text"
            )
            truncated = getattr(response, "stop_reason", None) == "max_tokens"
            if truncated and attempt == 1:
                logger.warning(
                    "Claude vision response hit max_tokens=%d on %s; retrying "
                    "once with max_tokens=%d",
                    budget, image_path.name, budget * 2,
                )
                continue
            parsed = _parse_complete_vision_response(
                response_text,
                backend="Claude vision",
                truncated=truncated,
                truncation_detail=(
                    f"stopped at max_tokens={budget} on {image_path.name} "
                    f"after {attempt} attempt(s); "
                    f"{len(expected_labels)} region label(s) expected"
                ),
            )
            break

        # One shared converter for both backends (#253). It was duplicated
        # here and in the local-VLM path, so the pixel-coordinate defect had
        # to be found and fixed twice.
        bbox_counts: Counter = Counter()

        def _norm_to_px(bbox_norm):
            px, disposition = _bbox_to_px(bbox_norm, w, h)
            if bbox_norm is not None:
                bbox_counts[disposition] += 1
            return px

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
        _log_bbox_dispositions(bbox_counts, image_path.name, self.name)
        return out


# ---------------------------------------------------------------------------
# Local VLM backend (Qwen2.5-VL)
# ---------------------------------------------------------------------------


_DEFAULT_LOCAL_VLM = DEFAULT_VISION_MODELS["vision-local"]

# Models smaller than 7B work on MPS / smaller GPUs.
_LOCAL_VLM_VARIANTS = {
    "Qwen/Qwen2.5-VL-7B-Instruct",
    "Qwen/Qwen2.5-VL-3B-Instruct",
    "Qwen/Qwen2.5-VL-2B-Instruct",
}


# The local model receives the same task as Claude, but phrased as a
# Qwen2.5-VL chat prompt.  Qwen's vision-language models respond to
# ``<|im_start|>system`` / ``<|im_start|>user`` templates automatically
# when run through the HuggingFace transformers chat pipeline.
_LOCAL_SYSTEM_PROMPT = _CLAUDE_SYSTEM_PROMPT  # identical task spec


class LocalVLMBackend(VisionBackend):
    """Pass 3b backend using a local Qwen2.5-VL model on CUDA / MPS / CPU.

    Loads the model and processor once at construction; all subsequent
    ``detect_figure_panels()`` calls run inference locally — no network
    required. On an H200, per-figure inference is ~1–3 s.

    Device selection follows the same ``CORPUS_DEVICE`` convention as
    ``pipeline.embeddings``. The default is auto-detect (cuda → mps → cpu).

    ``max_pixels`` controls the maximum image resolution fed to the
    model. Qwen2.5-VL resizes internally (min/max pixel budget), but
    for scientific figures the default 1280×28×28 grid is enough; going
    higher burns VRAM without improving label detection.

    ``max_new_tokens`` is a **floor**, not the budget (#253). The prompt
    asks for one six-field JSON object per panel, so the output length is
    a function of the panel count — which Pass 3b already knows from the
    caption before it calls the model. A fixed 1024 ran out mid-object on
    panel-rich figures, and a response that stops mid-object has no
    balanced ``{...}`` to extract, so the backend returned ``[]`` and the
    figure was recorded as *no-labels*: indistinguishable from a figure
    the model genuinely found nothing in. In one full reference run, ROI
    coverage fell from 47.8% at 2–3 panels to 13.6% at 10+, which is the
    signature of a fixed output budget rather than a vision failure. Bare
    plate discovery has no known target count, so it receives a conservative
    4096-token floor and still fails visibly if that cap is reached.
    """

    def __init__(
        self,
        model: str = _DEFAULT_LOCAL_VLM,
        *,
        device: Optional[str] = None,
        max_new_tokens: int = 1024,
        max_pixels: int = 1003520,  # 1280 * 28 * 28
        min_pixels: int = 3136,     # 4 * 28 * 28
    ):
        self._model_id = model
        self._max_new_tokens = max_new_tokens
        self._max_pixels = max_pixels
        self._min_pixels = min_pixels

        # Device selection — reuse embeddings.detect_device when
        # available; fall back to a simple torch probe.
        if device:
            self._device = device
        else:
            try:
                from .embeddings import detect_device
                self._device = detect_device()
            except ImportError:
                self._device = self._probe_device()

        logger.info("Loading local VLM %s on device=%s", model, self._device)

        try:
            from transformers import (
                Qwen2_5_VLForConditionalGeneration,
                AutoProcessor,
            )
        except ImportError as e:
            raise VisionBackendError(
                "transformers >= 4.45 is required for the local VLM backend "
                "(pip install transformers>=4.45 qwen-vl-utils torch "
                "accelerate)"
            ) from e

        try:
            import torch
            dtype = torch.bfloat16 if self._device == "cuda" else torch.float32
            self._processor = AutoProcessor.from_pretrained(
                model,
                min_pixels=self._min_pixels,
                max_pixels=self._max_pixels,
            )
            self._model = Qwen2_5_VLForConditionalGeneration.from_pretrained(
                model,
                torch_dtype=dtype,
                device_map=(self._device if self._device == "cuda" else None),
            )
            if self._device != "cuda":
                self._model = self._model.to(self._device)
        except Exception as e:
            raise VisionBackendError(
                f"Could not load local VLM {model!r} on device={self._device}: {e}"
            ) from e

    panel_mode = "vision-local"

    @cached_property
    def producer(self):
        result = vision_producer(self.panel_mode, self._model_id,
                                 resolved_revision=getattr(self._model.config, "_commit_hash", None))
        result["generation"] = {"max_new_tokens": self._max_new_tokens,
                                "max_pixels": self._max_pixels, "min_pixels": self._min_pixels}
        return result

    def _token_budget(self, expected_labels) -> int:
        """Output-token budget for one figure, scaled by its panel count.

        The configured ``max_new_tokens`` is the floor, so a caller that
        raised it keeps what they asked for and small figures are unaffected.
        A bare-plate discovery request has no caption-derived count to scale
        by, so the shared helper supplies its conservative discovery floor.
        """
        return _vision_token_budget(self._max_new_tokens, expected_labels)

    @staticmethod
    def _probe_device() -> str:
        try:
            import torch
        except ImportError:
            return "cpu"
        if torch.cuda.is_available():
            return "cuda"
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return "mps"
        return "cpu"

    @property
    def name(self) -> str:
        # Short identifier — just the model's basename.
        return f"vision:{self._model_id.split('/')[-1].lower()}"

    def detect_figure_panels(
        self,
        image_path: Path,
        caption_text: str,
        expected_labels: List[str],
    ) -> List[Dict]:
        import torch
        from PIL import Image

        try:
            img = Image.open(image_path).convert("RGB")
            w, h = img.size
        except Exception as e:
            raise VisionBackendError(
                f"could not read image {image_path}: {e}"
            ) from e

        user_text = _vision_user_text(caption_text, expected_labels, w, h)

        messages = [
            {"role": "system", "content": _LOCAL_SYSTEM_PROMPT},
            {"role": "user", "content": [
                {"type": "image", "image": img},
                {"type": "text", "text": user_text},
            ]},
        ]

        try:
            text_prompt = self._processor.apply_chat_template(
                messages, tokenize=False, add_generation_prompt=True,
            )
            from qwen_vl_utils import process_vision_info
            image_inputs, video_inputs = process_vision_info(messages)
            inputs = self._processor(
                text=[text_prompt],
                images=image_inputs,
                videos=video_inputs,
                padding=True,
                return_tensors="pt",
            ).to(self._model.device)

        except Exception as e:
            raise VisionBackendError(
                f"Local VLM input preparation failed on {image_path.name}: {e}"
            ) from e

        initial_budget = self._token_budget(expected_labels)
        for attempt, budget in enumerate(
            (initial_budget, initial_budget * 2), start=1,
        ):
            try:
                with torch.no_grad():
                    output_ids = self._model.generate(
                        **inputs,
                        max_new_tokens=budget,
                        do_sample=False,
                    )
                # Strip the prompt tokens to get only the generated response.
                generated = output_ids[:, inputs.input_ids.shape[1]:]
                response_text = self._processor.batch_decode(
                    generated, skip_special_tokens=True,
                )[0]
            except Exception as e:
                raise VisionBackendError(
                    f"Local VLM inference failed on {image_path.name}: {e}"
                ) from e

            # Qwen exposes no stop reason here. Reaching max_new_tokens is the
            # truncation signal. A response that happens to emit EOS on the
            # final slot gets one conservative retry; it can never be accepted
            # as a possibly partial answer (#269).
            truncated = int(generated.shape[1]) >= budget
            if truncated and attempt == 1:
                logger.warning(
                    "Local VLM response hit max_new_tokens=%d on %s; retrying "
                    "once with max_new_tokens=%d",
                    budget, image_path.name, budget * 2,
                )
                continue
            parsed = _parse_complete_vision_response(
                response_text,
                backend="Local VLM",
                truncated=truncated,
                truncation_detail=(
                    f"reached max_new_tokens={budget} on {image_path.name} "
                    f"after {attempt} attempt(s); "
                    f"{len(expected_labels)} region label(s) expected"
                ),
            )
            break

        # Convert normalized bboxes to pixel coordinates — same logic as
        # the Claude backend.
        # One shared converter for both backends (#253). It was duplicated
        # here and in the local-VLM path, so the pixel-coordinate defect had
        # to be found and fixed twice.
        bbox_counts: Counter = Counter()

        def _norm_to_px(bbox_norm):
            px, disposition = _bbox_to_px(bbox_norm, w, h)
            if bbox_norm is not None:
                bbox_counts[disposition] += 1
            return px

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
        _log_bbox_dispositions(bbox_counts, image_path.name, self.name)
        return out


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def get_vision_backend(backend: str = "claude", **kwargs) -> VisionBackend:
    """Construct a vision backend by name.

    ``"claude"`` — Anthropic Claude API (default for development).
    ``"local"``  — Qwen2.5-VL on CUDA / MPS / CPU (for production on
    Bouchet; no network or per-call cost).
    """
    backend = (backend or "claude").lower()
    if backend == "claude":
        return ClaudeVisionBackend(**kwargs)
    if backend == "local":
        return LocalVLMBackend(**kwargs)
    raise ValueError(f"Unknown vision backend: {backend!r}")

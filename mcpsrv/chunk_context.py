"""Bounded projections of immutable chunk context and build receipts (#319)."""
from __future__ import annotations

from itertools import islice
import json

MAX_CONTEXT_ROW_BYTES = 8192
MAX_OPTIONAL_CONTEXT_BYTES = 65536
MAX_CONTEXT_RECORDS = 6
MAX_NESTED_RECORDS = 6
MAX_FIELDS = 16
MAX_DEPTH = 4
MAX_NODES = 128
_CONTEXT_LISTS = ("source_items", "text_integrity", "tables", "key_branches")
_PROSE = {"original", "replacement", "before", "after", "text", "source_text",
          "original_lead", "original_destination", "heading", "caption",
          "original_native_token", "decoded_hint", "original_tokens"}
_TREATMENT_FIELDS = ("status", "name", "rank", "heading", "heading_ref", "heading_page",
                     "evidence", "name_evidence_ref")


def context_bytes(value):
    """Compact UTF-8 JSON bytes, excluding MCP transport/pretty-print overhead."""
    return len(json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))


def _preview(value, limit):
    return {"preview": value[:limit], "char_count": len(value),
            "preview_charspan": [0, min(limit, len(value))], "truncated": len(value) > limit}


def _counts(available, returned):
    return {"available": available, "returned": returned, "truncated": returned < available}


def _project(value, key, preview_limit, state, depth=0):
    """Bound every nested observation without scanning omitted list contents."""
    state["nodes"] += 1
    if state["nodes"] > MAX_NODES or depth > MAX_DEPTH:
        state["truncated"] = True
        return {"omitted": True, "reason": "projection_complexity_limit"}
    if isinstance(value, str):
        if key in _PROSE or len(value) > 256:
            result = _preview(value, preview_limit if key in _PROSE else 256)
            state["truncated"] |= result["truncated"]
            return result
        return value
    if value is None or isinstance(value, (bool, int, float)):
        # Malformed artifacts must not hide arbitrarily large integer strings.
        if isinstance(value, int) and value.bit_length() > 128:
            state["truncated"] = True
            return {"omitted": True, "reason": "oversized_number"}
        return value
    if isinstance(value, list):
        result = []
        for item in islice(value, MAX_NESTED_RECORDS):
            if state["nodes"] >= MAX_NODES:
                break
            result.append(_project(item, key, preview_limit, state, depth + 1))
        if len(result) < len(value):
            state["truncated"] = True
        return result
    if isinstance(value, dict):
        result = {}
        for field, item in islice(value.items(), MAX_FIELDS):
            if state["nodes"] >= MAX_NODES or len(field) > 64:
                state["truncated"] = True
                continue
            result[field] = _project(item, field, preview_limit, state, depth + 1)
            if isinstance(item, list):
                shown = len(result[field]) if isinstance(result[field], list) else 0
                if field not in {"charspan", "source_bbox", "base_bbox", "bbox"} or shown < len(item):
                    result[field + "_scope"] = _counts(len(item), shown)
        shown = sum(field in result for field in islice(value, MAX_FIELDS))
        if shown < len(value):
            result["_field_scope"] = _counts(len(value), shown)
            state["truncated"] = True
        return result
    state["truncated"] = True
    return {"omitted": True, "reason": "unsupported_context_value"}


class ContextProjection:
    """One response's optional evidence budget; never removes requested rows.

    Every row has its own complete byte ceiling for these newly added fields.
    The shared ceiling covers optional evidence arrays. Required status/count
    notices necessarily scale with the pre-existing response row count.
    """

    def __init__(self, *, with_text):
        self.preview_limit = 96 if with_text else 32
        self.optional_remaining = MAX_OPTIONAL_CONTEXT_BYTES

    def project(self, chunk):
        context = chunk.get("treatment_context")
        if not isinstance(context, dict) or not context:
            context = {"status": "unavailable", "name": None}
        state = {"nodes": 0, "truncated": False}
        treatment = {key: _project(context[key], key, self.preview_limit, state)
                     for key in _TREATMENT_FIELDS if key in context}
        treatment.setdefault("status", "unavailable")
        treatment.setdefault("name", None)
        section = _project(chunk.get("section_type"), "section_type", self.preview_limit, state)
        projection = {
            "policy": "bounded_context_previews_v1",
            "row_byte_limit": MAX_CONTEXT_ROW_BYTES,
            "optional_response_byte_limit": MAX_OPTIONAL_CONTEXT_BYTES,
            "preview_char_limit": self.preview_limit,
            "truncated": state["truncated"],
            "records": {},
        }
        result = {"treatment_context": treatment, "section_type": section,
                  "source_items": [], "context_projection": projection}
        sources = {}
        for field in _CONTEXT_LISTS:
            values = chunk.get(field)
            if not isinstance(values, list):
                values = []
            if field == "source_items" or field in chunk:
                sources[field] = values
                result[field] = []
                projection["records"][field] = _counts(len(values), 0)
        # Producer treatment names/refs are small. Fail visibly closed even for
        # malformed huge nested values in those fields, retaining row identity.
        if context_bytes(result) > MAX_CONTEXT_ROW_BYTES // 2:
            result["treatment_context"] = {"status": "unavailable", "name": None,
                                           "projection_omitted": True}
            result["section_type"] = None
            projection["truncated"] = True
        remaining = min(MAX_CONTEXT_ROW_BYTES - context_bytes(result), self.optional_remaining)
        used = 0
        for field, values in sources.items():
            for value in islice(values, MAX_CONTEXT_RECORDS):
                record_state = {"nodes": 0, "truncated": False}
                record = _project(value, field, self.preview_limit, record_state)
                cost = context_bytes(record) + (1 if result[field] else 2)  # comma or array brackets
                if cost > remaining:
                    projection["truncated"] = True
                    break
                result[field].append(record)
                remaining -= cost
                used += cost
                projection["records"][field] = _counts(len(values), len(result[field]))
                projection["truncated"] |= record_state["truncated"]
            if len(result[field]) < len(values):
                projection["truncated"] = True
        # Count/status digits can grow after records are added. Leave exact
        # per-row enforcement to this final bounded check, omitting whole
        # records rather than silently cutting a JSON object or source span.
        for field in reversed(sources):
            while result[field] and context_bytes(result) > MAX_CONTEXT_ROW_BYTES:
                removed = result[field].pop()
                used -= context_bytes(removed) + (1 if result[field] else 2)
                projection["records"][field] = _counts(len(sources[field]), len(result[field]))
                projection["truncated"] = True
        self.optional_remaining -= used
        return result

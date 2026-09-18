"""Conservative, auditable dispositions for extracted reference evidence (#313)."""
from __future__ import annotations

import json
import re

PRODUCER = "reference-quality-v2"


def create_schema(conn):
    conn.execute("""CREATE TABLE IF NOT EXISTS reference_observation_quality (
        observation_id TEXT PRIMARY KEY REFERENCES reference_observations(observation_id),
        disposition TEXT NOT NULL,
        reasons_json TEXT NOT NULL,
        producer_version TEXT NOT NULL
    )""")


def classify(conn, ref):
    """Quarantine recognizable caption/legend forms, never missing fields alone."""
    title = (ref.get("title") or "").strip()
    text = " ".join(str(ref.get(key) or "") for key in ("raw", "title", "journal"))
    # An explicit publication identifier/year is counter-evidence. This small
    # rule set does not attempt to adjudicate ambiguous bibliographic forms.
    if ref.get("doi") or ref.get("year"):
        return "usable", []
    # A curated library entry is independent evidence that this title really
    # is a publication. Its unusual typography alone cannot quarantine it.
    if title and conn.execute("SELECT 1 FROM works WHERE title=? AND bib_imported_at IS NOT NULL", (title,)).fetchone():
        return "usable", [{"code": "curated_title_counterevidence"}]
    reasons = []
    definitions = re.findall(r"\b[A-Za-z]+(?:[.-][A-Za-z]+)+\s*={1,2}\s*[A-Za-z]", text)
    if len(definitions) >= 2:
        reasons.append({"code": "multiple_legend_definitions"})
    panel_run = re.search(r"\bA\s*,\s*B\s*,\s*C\b", text)
    panel_views = panel_run and re.search(r"\b(?:views?|panels?)\b", text, re.IGNORECASE)
    figure_anchor = re.search(r"\b(?:fig(?:ure)?\.?|fic\.|plate)\s+(?:[IVXLCDM]+|\d+)\b", text, re.IGNORECASE)
    if panel_views and figure_anchor:
        reasons.append({"code": "caption_panel_view_description"})
    plate = re.compile(r"^(?:plate|fig(?:ure)?\.?|fic\.)\s+(?:[IVXLCDM]+|\d+)\b", re.IGNORECASE)
    raw = (ref.get("raw") or title or ref.get("journal") or "").strip()
    if plate.match(raw) and not ref.get("authors"):
        reasons.append({"code": "standalone_plate_or_figure_label"})
    if reasons:
        return "quarantined_fragment", reasons
    if panel_views:
        return "review_needed", [{"code": "possible_panel_description", "requires_source_review": True}]
    return "usable", []


def record(conn, observation_id, disposition, reasons):
    conn.execute("INSERT OR REPLACE INTO reference_observation_quality VALUES (?,?,?,?)",
                 (observation_id, disposition, json.dumps(reasons, sort_keys=True, ensure_ascii=False), PRODUCER))

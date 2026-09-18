"""Source-faithful table/key serialization and evidence-only word boundaries.

Build-plane helpers for #307, #308 and #334. A merged cell is one observation,
not a value to repeat across the grid. Missing spaces are restored only where
an independently extracted PDF line contains exactly the same letters with
spaces; source typography and unknown boundaries stay unchanged.
"""
from __future__ import annotations

from collections import defaultdict
from html import escape
import re

from docling_core.transforms.chunker.hierarchical_chunker import (
    ChunkingDocSerializer, ChunkingSerializerProvider,
)
from docling_core.transforms.serializer.base import BaseTableSerializer
from docling_core.transforms.serializer.common import create_ser_result
from docling_core.transforms.serializer.markdown import MarkdownDocSerializer
from docling_core.types.doc import TableItem
from docling_core.types.doc.common.meta import BaseMeta, FloatingMeta

TABLE_STRUCTURE_POLICY = "logical-cells-key-geometry-source-spaces-v1"
_LONG_RUN = re.compile(r"[^\W\d_]{20,}", re.UNICODE)
_LEADERS = re.compile(r"(?:\.\s*){3,}")


def _metadata(item, key, value):
    if item.meta is None:
        item.meta = FloatingMeta() if isinstance(item, TableItem) else BaseMeta()
    setattr(item.meta, key, value)


def _box(bbox, height):
    """Return top-left coordinates without changing stored provenance."""
    b = bbox.to_top_left_origin(page_height=height)
    return b.l, b.t, b.r, b.b


def restore_source_spaces(text, source_lines):
    """Insert spaces supported by an exact contiguous source-line match.

    Never splits a word by a lexicon, case, language, or estimated character gap.
    Ambiguous repeated matches must agree on every boundary. Newline boundaries
    are deliberately excluded: a source line break can be inside a compound.
    """
    evidence = []
    replacements = []
    for match in _LONG_RUN.finditer(text):
        token = match.group()
        options = []
        for source in source_lines:
            compact = "".join(c for c in source if not c.isspace())
            positions = [i for i, c in enumerate(source) if not c.isspace()]
            start = compact.find(token)
            while start >= 0:
                end = start + len(token)
                # Do not match the middle of a different source word.
                left = positions[start]
                right = positions[end - 1] + 1
                if (left == 0 or not source[left - 1].isalpha()) and (
                    right == len(source) or not source[right].isalpha()
                ):
                    candidate = " ".join(source[left:right].split())
                    options.append(candidate)
                start = compact.find(token, start + 1)
        if options and len(set(options)) == 1 and " " in options[0]:
            replacements.append((match.start(), match.end(), options[0]))
            evidence.append({"original": token, "replacement": options[0],
                             "route": "pdf_line_exact_letters"})
    for start, end, value in reversed(replacements):
        text = text[:start] + value + text[end:]
    return text, evidence


def logical_rows(table, doc=None, doc_serializer=None, **kwargs):
    """Return each actual cell once, with source offsets/spans intact."""
    rows = defaultdict(list)
    for index, cell in enumerate(table.data.table_cells):
        value = cell._get_text(doc=doc, **({"doc_serializer": doc_serializer}
                                         if doc_serializer is not None else {}), **kwargs)
        rows[cell.start_row_offset_idx].append({
            "cell_index": index, "row": cell.start_row_offset_idx,
            "column": cell.start_col_offset_idx,
            "row_span": cell.row_span, "col_span": cell.col_span,
            "column_header": cell.column_header, "row_header": cell.row_header,
            "text": value,
        })
    return [(row, sorted(cells, key=lambda c: (c["column"], c["cell_index"])))
            for row, cells in sorted(rows.items())]


def is_identification_key(rows):
    """Require a sequence of numbered narrative leads, not just numeric data."""
    numbered = []
    narrative = 0
    for _, cells in rows:
        text = " ".join(c["text"] for c in cells).strip()
        if len(re.findall(r"[^\W\d_]+", text)) >= 6:
            narrative += 1
            if m := re.match(r"^(\d{1,3})(?:\s|\.)", text):
                numbered.append(int(m[1]))
    return (len(rows) >= 4 and narrative >= len(rows) * .7
            and len(numbered) >= 2 and numbered == sorted(set(numbered)))


def serialized_rows(table, doc=None, doc_serializer=None, **kwargs):
    rows = logical_rows(table, doc, doc_serializer, **kwargs)
    key = is_identification_key(rows)
    result = []
    for row, cells in rows:
        if key:
            # Leaders carry alignment, not linguistic content. Their removal
            # prevents them alone consuming the token budget of a key branch.
            text = _LEADERS.sub(" … ", " ".join(c["text"] for c in cells))
            text = " ".join(text.split())
        else:
            columns = [""] * table.data.num_cols
            for cell in cells:
                value = " ".join(cell["text"].split()).replace("|", r"\|")
                col = cell["column"]
                if 0 <= col < len(columns):
                    columns[col] = (columns[col] + " " + value).strip()
            text = "| " + " | ".join(columns) + " |"
        if text.strip("| "):
            result.append((row, text, cells))
    return result, key


class SourceTableSerializer(BaseTableSerializer):
    """Emit logical rows, never DataFrame grid expansion or triplet prose."""

    def serialize(self, *, item, doc_serializer, doc, **kwargs):
        parts = []
        caption = doc_serializer.serialize_captions(item=item, **kwargs)
        if caption.text:
            parts.append(caption)
        if item.self_ref not in doc_serializer.get_excluded_refs(**kwargs):
            rows, _ = serialized_rows(item, doc, doc_serializer, **kwargs)
            if rows:
                parts.append(create_ser_result(
                    text="\n".join(text for _, text, _ in rows), span_source=item))
        return create_ser_result(text="\n\n".join(p.text for p in parts), span_source=parts)

    def get_header_and_body_lines(self, *, table_text, **kwargs):
        # Keep one logical row (including a whole key branch) per split unit.
        # Repeated header text would inflate annotation mention counts; its
        # source-cell context belongs in chunk metadata instead.
        return [], table_text.splitlines(keepends=True)


class SourceMarkdownTableSerializer(SourceTableSerializer):
    """HTML spans retain the actual table structure in the Markdown artifact."""

    def serialize(self, *, item, doc_serializer, doc, **kwargs):
        rows = logical_rows(item, doc, doc_serializer, **kwargs)
        if is_identification_key(rows):
            return super().serialize(item=item, doc_serializer=doc_serializer, doc=doc, **kwargs)
        parts = []
        caption = doc_serializer.serialize_captions(item=item, **kwargs)
        if caption.text:
            parts.append(caption)
        if item.self_ref not in doc_serializer.get_excluded_refs(**kwargs):
            covered = set()
            lines = ["<table>"]
            for row, cells in rows:
                values = []
                col = 0
                for cell in cells:
                    while col < cell["column"]:
                        if (row, col) not in covered:
                            values.append("<td></td>")
                        col += 1
                    tag = "th" if cell["column_header"] or cell["row_header"] else "td"
                    attributes = "".join(f' {key}="{cell[field]}"' for key, field in
                                         (("rowspan", "row_span"), ("colspan", "col_span"))
                                         if cell[field] > 1)
                    values.append(f"<{tag}{attributes}>{escape(cell['text'])}</{tag}>")
                    for r in range(row, row + cell["row_span"]):
                        for c in range(cell["column"], cell["column"] + cell["col_span"]):
                            covered.add((r, c))
                    col = cell["column"] + cell["col_span"]
                while col < item.data.num_cols:
                    if (row, col) not in covered:
                        values.append("<td></td>")
                    col += 1
                lines.append("<tr>" + "".join(values) + "</tr>")
            lines.append("</table>")
            parts.append(create_ser_result(text="\n".join(lines), span_source=item))
        return create_ser_result(text="\n\n".join(p.text for p in parts), span_source=parts)


def _block_build_metadata(serializer, document):
    # Docling otherwise serializes every custom metadata value as prose,
    # repeating raw repair histories and source labels in annotation text.
    for item in [*document.texts, *document.tables, *document.pictures, *document.groups]:
        if item.meta:
            serializer.params.blocked_meta_names.update(
                name for name in item.meta.model_dump() if name.startswith("corpus__"))
    return serializer


class SourceTableSerializerProvider(ChunkingSerializerProvider):
    def get_serializer(self, doc):
        return _block_build_metadata(ChunkingDocSerializer(
            doc=doc, table_serializer=SourceTableSerializer()), doc)


def export_source_markdown(document):
    """Use the same logical cells in text.json and chunks.json."""
    return _block_build_metadata(MarkdownDocSerializer(
        doc=document, table_serializer=SourceMarkdownTableSerializer()), document).serialize().text


def table_chunk_metadata(items, text, document):
    """Attach row/cell identity and continuation without repeating source text.

    Partial rows are explicit; a client must not treat their pieces as a whole
    branch. Source headers/couplet numbers live here rather than in annotation
    text. Source cell offsets and spans remain available in docling_doc.json.
    """
    tables = []
    for table in items:
        # HybridChunker's short-chunk copy can downcast TableItem to DocItem;
        # the source reference still resolves the complete structured table.
        if not isinstance(table, TableItem):
            try:
                table = table.get_ref().resolve(document)
            except (AttributeError, IndexError, KeyError, ValueError):
                continue
        if not isinstance(table, TableItem):
            continue
        rows, key = serialized_rows(table, document)
        selected = [(row, value, cells) for row, value, cells in rows if value in text]
        selected_ids = {row for row, _, _ in selected}
        caption_texts = [r.resolve(document).text for r in table.captions
                         if hasattr(r.resolve(document), "text")]
        partial = not selected or any(
            line.strip() and not any(line.strip() in value for _, value, _ in selected)
            and not any(line.strip() in caption for caption in caption_texts)
            for line in text.splitlines() if key or line.strip().startswith("|"))
        current = None
        couplets = {}
        for row, value, _ in rows:
            if m := re.match(r"^(\d{1,3})(?:\s|\.)", value):
                current = m[1]
            if row in selected_ids and current:
                couplets[str(row)] = current
        tables.append({
            "table_ref": table.self_ref, "kind": "identification_key" if key else "table",
            "row_indices": sorted(selected_ids), "total_rows": len(rows),
            "complete_rows": not partial, "continuation": len(selected) < len(rows),
            "source_cells": [{k: v for k, v in cell.items() if k != "text"}
                             for _, _, cells in selected for cell in cells],
            "column_headers": [cell for _, _, cells in rows for cell in cells
                               if cell["column_header"]],
            "spanning_context": [cell for row, _, cells in rows for cell in cells
                                 if row not in selected_ids and cell["row_span"] > 1
                                 and any(row <= r < row + cell["row_span"] for r in selected_ids)],
            "couplets": couplets if key else {},
        })
    return tables


def _attach_key_destinations(document):
    """Bind detached right-aligned endpoints to the unique lead on that line.

    Only explicit key sections and their continuation pages qualify. A nearby
    word alone is insufficient. Original observations remain in item metadata.
    """
    active_key = None
    scoped = []
    for item, _ in document.iterate_items():
        value = getattr(item, "text", "")
        if re.search(r"\bkey\s+to\b", value, re.I) and len(value) < 160:
            active_key = value
        elif str(item.label) == "section_header":
            active_key = None
        if active_key and value and getattr(item, "prov", None):
            scoped.append((item, active_key))
    pairs = []
    for endpoint, heading in scoped:
        target = endpoint.text.strip()
        if not re.fullmatch(r"[\w?.-]+(?:\s+[\w?.-]+){0,3}", target) or len(target) > 60:
            continue
        ep = endpoint.prov[0]
        page = document.pages[ep.page_no]
        ex0, ey0, _, ey1 = _box(ep.bbox, page.size.height)
        candidates = []
        for lead, lead_heading in scoped:
            if lead is endpoint or heading != lead_heading or not lead.prov:
                continue
            lp = lead.prov[0]
            if lp.page_no != ep.page_no or not re.search(r"\.\s*\.\s*[.\s]*$", lead.text):
                continue
            _, ly0, lx1, ly1 = _box(lp.bbox, page.size.height)
            if (ex0 >= lx1 - 1 and ex0 - lx1 < page.size.width * .4
                    and abs(ey1 - ly1) <= max(3, (ey1 - ey0) * .4)
                    and ey0 >= ly0 - 2):
                candidates.append(lead)
        if len(candidates) == 1:
            pairs.append((candidates[0], endpoint, heading))
    # Multiple endpoints competing for one lead are unresolved, not guessed.
    counts = defaultdict(int)
    for lead, _, _ in pairs:
        counts[lead.self_ref] += 1
    report = []
    for lead, endpoint, heading in pairs:
        if counts[lead.self_ref] != 1:
            continue
        observation = {"heading": heading, "lead_ref": lead.self_ref,
                       "destination_ref": endpoint.self_ref,
                       "original_lead": lead.text, "original_destination": endpoint.text,
                       "lead_provenance": [p.model_dump(mode="json") for p in lead.prov],
                       "destination_provenance": [p.model_dump(mode="json") for p in endpoint.prov],
                       "association": "unique_right_aligned_final_line",
                       "status": "geometry_verified_spelling_unverified"}
        lead.text = lead.text.rstrip() + " " + endpoint.text.strip()
        lead.prov = list(lead.prov) + list(endpoint.prov)
        _metadata(lead, "corpus__key_branch", observation)
        _metadata(endpoint, "corpus__key_destination", {"absorbed_by": lead.self_ref})
        # Keep the observation node and its original text/identity. Detaching
        # an orphan from the Docling tree can renumber unrelated references on
        # reload; an empty consumed node emits no second source occurrence.
        endpoint.text = ""
        report.append(observation)
    return report


def prepare_table_structure(document, pdf_path):
    """Prepare build observations; return a serializable, reviewable receipt."""
    import fitz

    from .source_spaces import (corroborate_source_gaps, geometric_space_candidate,
                                source_spacing_producer)

    producer = source_spacing_producer()
    report = {"policy": TABLE_STRUCTURE_POLICY, "space_repairs": [],
              "unresolved_long_runs": [], "key_associations": [],
              "source_spacing_producer": producer, "source_gap_observations": []}
    with fitz.open(str(pdf_path)) as source:
        # A page-at-a-time traversal avoids retaining a monograph's text cells.
        objects = defaultdict(list)
        for item in document.texts:
            if item.prov and _LONG_RUN.search(item.text):
                objects[item.prov[0].page_no].append((item, item, item.prov[0].bbox, None))
        for table in document.tables:
            if table.prov:
                for index, cell in enumerate(table.data.table_cells):
                    if cell.bbox and _LONG_RUN.search(cell.text):
                        objects[table.prov[0].page_no].append((table, cell, cell.bbox, index))
        for page_no, items in sorted(objects.items()):
            if not 1 <= page_no <= len(source):
                continue
            page = source[page_no - 1]
            lines = [(line["bbox"], "".join(s["text"] for s in line["spans"]))
                     for block in page.get_text("dict", flags=fitz.TEXTFLAGS_DICT & ~fitz.TEXT_PRESERVE_IMAGES)["blocks"]
                     for line in block.get("lines", [])]
            candidates = []
            for block in page.get_text("rawdict", flags=fitz.TEXTFLAGS_RAWDICT & ~fitz.TEXT_PRESERVE_IMAGES)["blocks"]:
                for line in block.get("lines", []):
                    candidate = geometric_space_candidate(line)
                    if candidate and any(choice["original"] in item.text
                                         for choice in candidate["choices"]
                                         for _, item, _, _ in items):
                        candidates.append(candidate)
            verified, observations = corroborate_source_gaps(page, candidates, producer)
            report["source_gap_observations"].extend(
                {**entry, "page": page_no} for entry in observations)
            confirmed = {choice["original"] for entry in observations
                         for choice in entry.get("accepted", [])}
            for value, entry in zip(verified, [e for e in observations if e.get("accepted")]):
                lines.append((entry["bbox"], value))
                if entry.get("isolated_font_runs"):
                    # Both OCR modes corroborate the internal gaps of these
                    # exact native runs. A captured item must already contain
                    # the complete run as a word: restore_source_spaces never
                    # adds a boundary at the surrounding font transition.
                    lines.extend((entry["bbox"], choice["candidate"])
                                 for choice in entry["accepted"])
            # A raw no-space native line does not compete with its verified
            # geometry+OCR representation when exact-letter matching below.
            lines = [(box, value) for box, value in lines
                     if value in verified or not any(token in value for token in confirmed)]
            for owner, item, bbox, index in items:
                x0, y0, x1, y1 = _box(bbox, page.rect.height)
                selected = [value for box, value in lines
                            if min(x1, box[2]) > max(x0, box[0])
                            and min(y1, box[3]) > max(y0, box[1])]
                original = item.text
                repaired, evidence = restore_source_spaces(original, selected)
                if evidence:
                    for proof in evidence:
                        gap_evidence = [
                            {"page": page_no, "crop_sha256": entry["crop_sha256"],
                             "source_charspan": [choice["start"], choice["end"]],
                             "original": choice["original"], "replacement": choice["candidate"]}
                            for entry in observations
                            if min(x1, entry["bbox"][2]) > max(x0, entry["bbox"][0])
                            and min(y1, entry["bbox"][3]) > max(y0, entry["bbox"][1])
                            for choice in entry.get("accepted", [])
                            if choice["original"] in proof["original"]
                            and choice["candidate"] in proof["replacement"]]
                        if gap_evidence:
                            proof["route"] = producer["policy"]
                            proof["source_gap_evidence"] = gap_evidence
                    item.text = repaired
                    observation = {"item_ref": owner.self_ref, "cell_index": index,
                                   "page": page_no, "bbox": bbox.model_dump(mode="json"),
                                   "original": original, "replacement": repaired,
                                   "evidence": evidence}
                    report["space_repairs"].append(observation)
                    history = list(getattr(owner.meta, "corpus__source_spaces", []) or [])
                    _metadata(owner, "corpus__source_spaces", history + [observation])
                remaining = [m.group() for m in _LONG_RUN.finditer(repaired)]
                if remaining:
                    report["unresolved_long_runs"].append({
                        "item_ref": owner.self_ref, "cell_index": index, "page": page_no,
                        "runs": remaining, "status": "no_source_space_evidence"})
    report["key_associations"] = _attach_key_destinations(document)
    return report

"""Materialized text-key branch coverage and chunk navigation (#308).

Source association belongs to the geometry producer. These helpers describe
how its already-associated text survived chunking; they never infer a name or
bind an unassociated destination from prose.
"""
from collections import defaultdict
import re

KEY_BRANCH_CONTEXT_POLICY = "source-key-branch-fragments-v1"


def key_branch_context(doc_items, chunk_text=None):
    records = []
    compact_chunk = " ".join(chunk_text.split()) if chunk_text is not None else None
    for item in doc_items:
        branch = getattr(getattr(item, "meta", None), "corpus__key_branch", None)
        if not branch:
            continue
        scope = {"coverage": "unknown", "complete_branch": None,
                 "continuation": None, "destination_in_text": None}
        record = {"item_ref": item.self_ref, "chunk_scope": scope}
        if compact_chunk is not None:
            target = " ".join(branch["original_destination"].split())
            # HybridChunker may downcast item metadata to DocItem, which has
            # provenance but no text. The geometry producer retains the exact
            # two source strings it joined, independently of that downcast.
            full_branch = " ".join((branch["original_lead"] + " " + target).split())
            complete = bool(full_branch) and bool(re.search(
                r"(?<!\w)" + re.escape(full_branch) + r"(?!\w)", compact_chunk))
            # A split single-item fragment can acquire a Markdown list prefix.
            # Require its exact remaining source text. Unalignable/mixed
            # excerpts remain unknown rather than estimating coverage.
            fragment = re.sub(r"^(?:[-*]|\d+[.)])\s+", "", compact_chunk)
            partial = bool(fragment) and any(c.isalpha() for c in fragment) and bool(re.search(
                r"(?<!\w)" + re.escape(fragment) + r"(?!\w)", full_branch))
            # Scope precedes verbose observations so a bounded projection
            # retains these compact semantics before source-text previews.
            scope["destination_in_text"] = bool(target and re.search(
                r"(?<!\w)" + re.escape(target) + r"(?!\w)", compact_chunk))
            if complete or partial:
                scope.update(coverage="complete" if complete else "partial",
                             complete_branch=complete, continuation=not complete)
        records.append({**record, **branch})
    return records


def link_key_fragments(chunks):
    """Link stored fragments by source item, without copying complete prose.

Only relationships already evidenced by the producer receive links. Navigation
is bounded to adjacent chunks even if one long branch spans many chunks.
"""
    branches = defaultdict(list)
    for chunk in chunks:
        for branch in chunk.get("key_branches", []):
            if "chunk_scope" in branch:
                branches[branch["item_ref"]].append((chunk["chunk_id"], branch))
    for fragments in branches.values():
        for index, (_, branch) in enumerate(fragments):
            branch["chunk_scope"].update({
                "fragment_index": index,
                "fragment_count": len(fragments),
                "previous_chunk_id": fragments[index - 1][0] if index else None,
                "next_chunk_id": fragments[index + 1][0] if index + 1 < len(fragments) else None,
            })

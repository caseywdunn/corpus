"""Materialize evidenced taxonomic treatment context, separately from mentions.

A taxon appearing in prose is not a treatment boundary. Only an explicit
heading (or a standalone new-species heading missed by layout detection) can
start one. Missing/ambiguous headings remain unknown.
"""
from __future__ import annotations

from copy import deepcopy
from difflib import SequenceMatcher
import re
import shutil
import subprocess

from .config import classify_section
from .source_layout import _label, item_bounds

TREATMENT_CONTEXT_POLICY = "source_treatments_v2"
_NAME = re.compile(r"^([A-Z][a-z]{2,})\s+([a-z][a-z/-]{2,})\b(.*)$", re.S)
_NEW_SPECIES = re.compile(r"(?:sp\s*\.?\s*nov\s*\.?|n\s*\.?\s*sp\s*\.?)\s*$", re.I)
_SECTION = re.compile(
    r"^(Diagnosis|Description(?:\s+of\s+(?:the\s+)?(?:holotype|paratypes?))?|"
    r"Material\s+examined|Holotype|Paratypes?|General\s+appearance|Etymology|Distribution|"
    r"Type\s+locality|Remarks|Notes|Nectosome|Siphosome|Pneumatophore)\s*(?:[.:]|$)", re.I)
_NON_TREATMENT = {"abstract", "introduction", "methods", "results", "discussion", "conclusion",
                  "references", "acknowledgements", "appendix"}


def _authority_names(text):
    """Require author-shaped tokens, not arbitrary prose before a year.

    This checks printed heading syntax only. It neither resolves an authority
    nor uses taxonomy to infer an omitted treatment name.
    """
    particles = {"de", "del", "della", "di", "da", "dos", "du", "van", "von", "der", "den", "le", "la"}
    text = re.sub(r"\s+et\s+al\.?$", "", text.strip())
    groups = re.split(r"\s+(?:and|&)\s+|\s*,\s*", text)
    for group in groups:
        tokens = group.split()
        if not tokens or not any(token[0].isupper() for token in tokens):
            return False
        for token in tokens:
            if token in particles:
                continue
            core = token.removesuffix(".")
            if not core or not core[0].isupper() or not all(c.isalpha() or c in "-'’" for c in core):
                return False
    return True


def _treatment_suffix(suffix):
    # Parenthesized original authorities are common. An unmatched parenthesis
    # or material after the authority remains unconfirmed.
    if suffix.startswith("(") and suffix.endswith(")"):
        suffix = suffix[1:-1].strip()
    new_species = _NEW_SPECIES.search(suffix)
    if new_species:
        authors = suffix[:new_species.start()].strip().removesuffix(",").strip()
        return not authors or _authority_names(authors)
    authority = re.fullmatch(r"(.+?)(?:,\s*|\s+)(?:17|18|19|20)\d{2}[a-z]?\.?", suffix)
    return bool(authority and _authority_names(authority[1]))


def recover_section_headings(document, pdf_path):
    """Re-read a narrowly suspect diagnostic heading from the rendered source.

    Similar spelling only selects a crop to inspect. It is never sufficient to
    replace text: independent rendered OCR must return the exact section word.
    """
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta
    report = {"method": "rendered_section_heading_v1", "repairs": [], "unresolved": []}
    executable = shutil.which("tesseract")
    with fitz.open(pdf_path) as pdf:
        for item in document.texts:
            if _label(item) != "section_header" or classify_section([item.text]):
                continue
            letters = re.sub(r"[^a-z]", "", item.text.lower())
            if not 6 <= len(letters) <= 14:
                continue
            options = sorted(((SequenceMatcher(None, letters, word).ratio(), word)
                              for word in ("diagnosis", "description")), reverse=True)
            if options[0][0] < .65 or options[0][0]-options[1][0] < .1:
                continue
            box = item_bounds(item, document)
            if not box or not 1 <= box[0] <= len(pdf):
                continue
            original = item.text
            observation = {"item_ref": item.self_ref, "page": box[0], "before": original,
                           "source_bbox": list(box[1:]), "method": report["method"]}
            confirmed = False
            if executable:
                rect = fitz.Rect(box[1:])+(-2,-2,2,2)
                pix = pdf[box[0]-1].get_pixmap(clip=rect, dpi=450)
                try:
                    result = subprocess.run([executable,"stdin","stdout","-l","eng","--psm","7"],
                                            input=pix.tobytes("png"), capture_output=True, timeout=10)
                    reading = re.sub(r"[^a-z]", "", result.stdout.decode("utf-8", errors="replace").lower())
                    confirmed = result.returncode == 0 and reading == options[0][1]
                except (OSError, subprocess.TimeoutExpired):
                    pass
            if confirmed:
                item.text = options[0][1].capitalize()
                observation.update({"after": item.text, "status": "repaired"})
                report["repairs"].append(observation)
            else:
                observation.update({"reason": "rendered_heading_not_confirmed", "status": "unresolved"})
                report["unresolved"].append(observation)
            if item.meta is None:
                item.meta = BaseMeta()
            item.meta.corpus__section_heading = observation
    return report


def _name_heading(item):
    text = " ".join(getattr(item, "text", "").split())
    match = _NAME.match(text)
    if not match or len(text) > 180:
        return None
    suffix = match[3].strip()
    if _label(item) != "section_header" and not _NEW_SPECIES.fullmatch(suffix):
        return None
    if not _treatment_suffix(suffix):
        return None
    return match[1]+" "+match[2], text, suffix


def _corroborated_name(candidate, suffix, following):
    """Resolve a damaged epithet only from the same nearby name and authority."""
    if "/" not in candidate:
        return candidate, None
    genus = candidate.split()[0]
    year = re.search(r"\b(?:17|18|19|20)\d{2}\b", suffix)
    matches = []
    for item in following[:5]:
        if _label(item) in {"caption", "picture", "table", "page_header", "page_footer"}:
            continue
        text = " ".join(getattr(item, "text", "").split())
        match = _NAME.match(text)
        if not match or match[1] != genus or not year or year[0] not in match[3][:80]:
            continue
        name = match[1]+" "+match[2]
        changes = [op for op in SequenceMatcher(None, candidate, name, autojunk=False).get_opcodes()
                   if op[0] != "equal"]
        if len(changes) == 1:
            _, a,b,c,d = changes[0]
            if candidate[a:b] == "/" and name[c:d] == "l":
                matches.append((name, item.self_ref))
    unique = {name for name, _ in matches}
    return matches[0] if len(unique) == 1 else (None, None)


def materialize_treatment_context(document):
    """Return item-reference metadata; do not alter literal taxon annotations."""
    items = [item for item, _ in document.iterate_items()]
    contexts = {}
    current = {"status": "unknown", "name": None}
    section_type = None
    section_evidence = None
    for i,item in enumerate(items):
        label = _label(item)
        text = " ".join(getattr(item, "text", "").split())
        box = item_bounds(item, document)
        if label in {"caption", "picture", "table", "page_header", "page_footer", "footnote"}:
            contexts[item.self_ref] = {"treatment_context": {"status": "unknown", "name": None},
                                       "section_type": None, "role": label}
            continue
        named = _name_heading(item)
        if named:
            name, heading, suffix = named
            name, corroboration = _corroborated_name(name, suffix, items[i+1:])
            current = {"status": "resolved" if name else "unknown", "name": name,
                       "heading": heading, "heading_ref": item.self_ref,
                       "heading_page": box[0] if box else None,
                       "evidence": "standalone_treatment_heading"}
            if corroboration:
                current["name_evidence_ref"] = corroboration
                current["evidence"] = "heading_and_repeated_name_authority"
            section_type = None
            section_evidence = None
        elif label == "section_header":
            cls = classify_section([text])
            if cls in _NON_TREATMENT or re.match(r"^(?:Genus|Family|Order|Subfamily|Suborder)\b", text):
                current = {"status": "unknown", "name": None}
            section_type = None
            section_evidence = None
        section = _SECTION.match(text)
        if section and (label == "section_header" or section.end() < len(text)):
            section_type = re.sub(r"\s+", "_", section[1].lower())
            section_evidence = getattr(getattr(item, "meta", None), "corpus__section_heading", None)
        contexts[item.self_ref] = {"treatment_context": deepcopy(current),
                                   "section_type": section_type, "role": label,
                                   "section_evidence": deepcopy(section_evidence)}
    # Captions may be reached through a picture serializer without appearing in
    # the body walk. Their names describe the illustration, not a treatment.
    for item in document.texts:
        if _label(item) == "caption":
            contexts[item.self_ref] = {"treatment_context": {"status": "unknown", "name": None},
                                       "section_type": None, "role": "caption"}
    return contexts


def chunk_source_context(doc_items, context_by_ref):
    """Collapse only unanimous treatment/section evidence across a chunk."""
    evidence = [context_by_ref.get(item.self_ref, {}) for item in doc_items]
    treatments = [e.get("treatment_context", {"status": "unknown", "name": None}) for e in evidence]
    treatment = treatments[0] if treatments and all(t == treatments[0] for t in treatments) else {"status": "unknown", "name": None}
    sections = {e.get("section_type") for e in evidence}
    section = next(iter(sections)) if len(sections) == 1 else None
    source_items = []
    text_integrity = []
    for entry in evidence:
        if entry.get("section_evidence") and entry["section_evidence"] not in text_integrity:
            text_integrity.append(entry["section_evidence"])
    key_branches = []
    for item in doc_items:
        meta = getattr(item, "meta", None)
        text_integrity.extend(getattr(meta, "corpus__scientific_text", []) or [])
        text_integrity.extend(getattr(meta, "corpus__text_encoding", []) or [])
        text_integrity.extend(getattr(meta, "corpus__native_text_recovery", []) or [])
        text_integrity.extend(getattr(meta, "corpus__surname_recovery", []) or [])
        heading_repair = getattr(meta, "corpus__section_heading", None)
        if heading_repair:
            text_integrity.append(heading_repair)
        branch = getattr(meta, "corpus__key_branch", None)
        if branch:
            key_branches.append({"item_ref": item.self_ref, **branch})
        for prov in getattr(item, "prov", []):
            source_items.append({"item_ref": item.self_ref, "page": prov.page_no,
                                 "bbox": prov.bbox.model_dump(mode="json"),
                                 "charspan": list(prov.charspan)})
    result = {"treatment_context": deepcopy(treatment), "section_type": section,
              "source_items": source_items}
    if text_integrity:
        result["text_integrity"] = text_integrity
    if key_branches:
        result["key_branches"] = key_branches
    return result


def context_merge_key(chunk, context_by_ref):
    """Never merge prose with captions or across an evidenced treatment boundary."""
    data = chunk_source_context(chunk.meta.doc_items, context_by_ref)
    treatment = data["treatment_context"]
    roles = frozenset(context_by_ref.get(item.self_ref, {}).get("role")
                      for item in chunk.meta.doc_items)
    media = bool(roles & {"caption", "picture", "table"})
    return (treatment.get("heading_ref"), data["section_type"],
            tuple(item.self_ref for item in chunk.meta.doc_items) if media else None)

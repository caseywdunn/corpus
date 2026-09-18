"""Source-backed citation groups and conservative author/year links (#309, #317).

Grobid expands grouped markers before serializing TEI. Its ref text is an
observation, not necessarily the words printed on the page. Recovering a
source span uses PDF coordinates; no repeated-author deletion rule is used.
"""
from __future__ import annotations

import re
from collections import OrderedDict
import unicodedata
from pathlib import Path

NS = {"tei": "http://www.tei-c.org/ns/1.0"}
_YEAR = re.compile(r"(?<!\d)((?:1[5-9]|20)\d{2})([a-z]?)(?!\d)")
_SEPARATOR = re.compile(r"^(?:[\s,;:()\[\]&–—-]|\b(?:and|or|[a-z])\b)*$")
_SUFFIXES = re.compile(r"\s*(?:([-–])\s*([a-z])\b|(?:,|and|&)\s*([a-z])\b)")


def _norm(text):
    return "".join(c for c in unicodedata.normalize("NFKD", text.casefold())
                   if c.isalnum())


class PdfCitationSource:
    """Read prepared-PDF characters once per cited page, without changing it."""

    def __init__(self, pdf_path: Path):
        import pymupdf
        self.doc = pymupdf.open(pdf_path)
        self.cache = OrderedDict()

    def close(self):
        self.doc.close()

    def page_chars(self, page):
        if page not in self.cache:
            chars = []
            import pymupdf
            flags = pymupdf.TEXTFLAGS_RAWDICT & ~pymupdf.TEXT_PRESERVE_IMAGES
            for block in self.doc[page - 1].get_text("rawdict", flags=flags)["blocks"]:
                for line in block.get("lines", []):
                    for span in line.get("spans", []):
                        chars.extend((c["c"], c["bbox"]) for c in span["chars"])
                    chars.append(("\n", None))
            self.cache[page] = chars
            if len(self.cache) > 4:
                self.cache.popitem(last=False)
        self.cache.move_to_end(page)
        return self.cache[page]

    def group_text(self, refs, raw_text=None):
        """Take the contiguous source interval between the group's ref boxes.

        An interval includes punctuation between boxes (e.g. a comma at a
        line end which Grobid omits from both refs). Per-page intervals avoid
        introducing running headers when a group straddles a page boundary.
        """
        page_boxes = {}
        try:
            for ref in refs:
                coords = ref.get("coords")
                if not coords:
                    return None
                for box in coords.split(";"):
                    page, x, y, w, h = map(float, box.split(","))
                    if page != int(page) or page < 1 or w <= 0 or h <= 0:
                        return None
                    page_boxes.setdefault(int(page), []).append((x, y, x + w, y + h))
            parts = []
            for page, boxes in page_boxes.items():
                chars = self.page_chars(page)
                hits = []
                for i, (_, bbox) in enumerate(chars):
                    if bbox is None:
                        continue
                    cy = (bbox[1] + bbox[3]) / 2
                    if any(min(bbox[2], x1) - max(bbox[0], x0) >= .25 * (bbox[2] - bbox[0])
                           and y0 - .1 <= cy <= y1 + .1
                           for x0, y0, x1, y1 in boxes):
                        hits.append(i)
                if not hits:
                    return None
                parts.append("".join(c for c, _ in chars[min(hits):max(hits) + 1]))
            source = " ".join(" ".join(parts).split())
        except (ValueError, IndexError, RuntimeError):
            return None
        raw = raw_text if raw_text is not None else " ".join("".join(r.itertext()) for r in refs)
        # Coordinate/text disagreement is reviewable, never an invitation to
        # substitute unrelated words. Whitespace and expansion differences
        # are expected, but the source and TEI must retain the same years.
        source_years = {m[1] for m in _YEAR.finditer(source)}
        raw_years = {m[1] for m in _YEAR.finditer(raw)}
        if not source or source_years != raw_years:
            return None
        source_key, raw_key = _norm(source), _norm(raw)
        cursor = iter(raw_key)
        if not all(any(c == candidate for candidate in cursor) for c in source_key):
            return None
        return source


def reference_evidence(root):
    """Retain surname lists and date suffixes which integer years discard."""
    refs = []
    for node in root.findall(".//tei:listBibl/tei:biblStruct", NS):
        xml_id = node.get("{http://www.w3.org/XML/1998/namespace}id")
        authors = node.findall("tei:analytic/tei:author", NS)
        if not authors:
            authors = node.findall("tei:monogr/tei:author", NS)
        surnames = [" ".join("".join(a.itertext()).split())
                    for author in authors for a in author.findall(".//tei:surname", NS)]
        date = node.find(".//tei:date[@type='published']", NS)
        date_text = "" if date is None else (date.get("when") or "".join(date.itertext()))
        year = _YEAR.search(date_text)
        raw = node.find("tei:note[@type='raw_reference']", NS)
        raw_text = "" if raw is None else "".join(raw.itertext())
        raw_year = _YEAR.search(raw_text)
        if year is None or not xml_id:
            continue
        suffix = year[2]
        if raw_year is not None and raw_year[1] == year[1]:
            suffix = raw_year[2] or suffix
        refs.append({"target": "#" + xml_id, "year": year[1], "suffix": suffix,
                     "surnames": surnames})
    return refs


def _author_key(text):
    text = re.sub(r"^\s*(?:and|see|e\.g\.)\b", "", text, flags=re.I)
    return _norm(re.sub(r"\b(?:and|et|al)\b", "", text, flags=re.I))


def _authors_match(author_text, reference):
    names = reference["surnames"]
    if not names:
        return False
    if re.search(r"\bet\s+al\b", author_text, re.I):
        return len(names) > 1 and _author_key(author_text) == _norm(names[0])
    return _author_key(author_text) == "".join(_norm(s) for s in names)


def _year_items(text):
    """Yield explicit years, including letter suffix lists/ranges, with spans."""
    consumed = 0
    for match in _YEAR.finditer(text):
        if match.start() < consumed:
            continue
        end = match.end()
        suffixes = [match[2]]
        if match[2]:
            while tail := _SUFFIXES.match(text, end):
                last = tail[2] or tail[3]
                if tail[1]:
                    if not suffixes[-1] <= last or ord(last) - ord(suffixes[-1]) > 10:
                        break
                    suffixes.extend(chr(c) for c in range(ord(suffixes[-1]) + 1, ord(last) + 1))
                else:
                    suffixes.append(last)
                end = tail.end()
        yield match.start(), end, match[1], list(dict.fromkeys(suffixes))
        consumed = end


def resolve_group(text, references):
    """Return complete source author/year spans, without guessing among ties."""
    if re.search(r"\d{4}\s*[-–]\s*\d{4}", text):
        return None  # A date range can describe one multi-volume work.
    records = []
    previous_end = 0
    author_start = author_end = None
    for start, end, year, suffixes in _year_items(text):
        between = text[previous_end:start]
        left = len(between) - len(between.lstrip(" \t\n,;:()[]&–—-"))
        right = len(between.rstrip(" \t\n,;:()[]&–—-"))
        explicit = between[left:right]
        if explicit:
            author_start, author_end = previous_end + left, previous_end + right
        if author_start is None:
            return None  # Numeric styles and author-outside-ref need their TEI observations.
        author = text[author_start:author_end]
        for suffix in suffixes:
            candidates = [r for r in references if r["year"] == year
                          and _authors_match(author, r)
                          and (not suffix or r["suffix"] == suffix)]
            if suffix and any(r["suffix"] == suffix for r in candidates):
                candidates = [r for r in candidates if r["suffix"] == suffix]
            target = candidates[0]["target"] if len(candidates) == 1 else None
            records.append({"target_xml_id": target, "surface": text[author_start:end],
                            "candidate_target_xml_ids": [r["target"] for r in candidates],
                            "citation_year": year + suffix,
                            "author_span": [author_start, author_end],
                            "year_span": [start, end],
                            "validation_status": "validated_author_year" if target else
                            ("ambiguous_author_year" if candidates else "unresolved_author_year")})
        previous_end = end
    return records or None


def _pieces(node):
    if node.text:
        yield node.text
    for child in node:
        if child.tag == "{%s}ref" % NS["tei"] and child.get("type") == "bibr":
            yield child
        else:
            yield from _pieces(child)
        if child.tail:
            yield child.tail


def paragraph_citations(para, references, source=None):
    """Render a paragraph with original source groups and span-linked records."""
    pieces = list(_pieces(para))
    rendered = []
    records = []
    i = 0
    while i < len(pieces):
        if isinstance(pieces[i], str):
            rendered.append(pieces[i])
            i += 1
            continue
        refs = [pieces[i]]
        raw_parts = ["".join(pieces[i].itertext())]
        j = i + 1
        while j < len(pieces):
            if isinstance(pieces[j], str):
                if (not _SEPARATOR.fullmatch(pieces[j]) or j + 1 >= len(pieces)
                        or isinstance(pieces[j + 1], str)):
                    break
                raw_parts.append(pieces[j])
                j += 1
            refs.append(pieces[j])
            raw_parts.append("".join(pieces[j].itertext()))
            j += 1
        raw = "".join(raw_parts)
        observed = [{"surface": "".join(r.itertext()).strip(),
                     "target_xml_id": r.get("target") or None,
                     "coords": r.get("coords")} for r in refs]
        recovered = source.group_text(refs, raw) if source is not None else None
        group = recovered if recovered is not None else raw
        # The whitespace boundary before a marker belongs to surrounding
        # prose; preserve it while normalizing the eventual stored paragraph.
        prefix = "".join(rendered)
        normalized_prefix = " ".join(prefix.split())
        group_start = len(normalized_prefix) + (1 if normalized_prefix and prefix[-1].isspace() else 0)
        group = " ".join(group.split())
        resolved = resolve_group(group, references)
        if resolved is not None and references:
            for rec in resolved:
                rec["tei_observations"] = observed
                rec["text_source"] = "pdf_coordinates" if recovered is not None else "tei"
                for key in ("author_span", "year_span"):
                    rec[key] = [n + group_start for n in rec[key]]
                records.append(rec)
        else:
            for observation in observed:
                if observation["surface"]:
                    records.append({"target_xml_id": observation["target_xml_id"],
                                    "surface": observation["surface"],
                                    "tei_observations": [observation],
                                    "text_source": "pdf_coordinates" if recovered is not None else "tei",
                                    "validation_status": "unverified_tei"})
        rendered.append(group)
        i = j
    return " ".join("".join(rendered).split()), records

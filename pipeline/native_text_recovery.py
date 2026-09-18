"""Retain damaged original-layer evidence across PDF preparation (#312).

Reversible decoding selects suspect tokens; it never supplies replacement text.
Two rendered regional OCR readings must agree, including the decoded accent.
Their source geometry then locates the corresponding prepared-PDF token before
any change to structured extraction. The library PDF is never modified.
"""
from __future__ import annotations

from collections import defaultdict
from difflib import SequenceMatcher
from functools import lru_cache
import hashlib
import importlib.metadata
import os
from pathlib import Path
import re
import shutil
import subprocess
import unicodedata

NATIVE_TEXT_RECOVERY_POLICY = "original-mojibake-regions-dual-raster-ocr-v1"
_DAMAGE = re.compile(r"(?:Ã[\x80-\xbf]|Â[\x80-\xbf])")
_DPIS = (300, 600)
_MAX_REGIONS_PER_PAGE = 16
_MAX_REGIONS_PER_DOCUMENT = 64
_TIMEOUT = 15


@lru_cache(maxsize=4)
def _installation(executable, modified, tessdata_prefix):
    del modified, tessdata_prefix
    try:
        version = subprocess.run([executable,"--version"],capture_output=True,text=True,
                                 check=True,timeout=5).stdout.strip()
        listing = subprocess.run([executable,"--list-langs"],capture_output=True,text=True,
                                 check=True,timeout=5).stdout
        location = re.search(r'"([^"]+)"', listing)
        directory = Path(location[1]) if location else None
        return version, directory
    except (OSError, subprocess.SubprocessError):
        return None, None


@lru_cache(maxsize=128)
def _model_hash(path, modified, size):
    del modified,size
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def native_text_recovery_producer(languages=None):
    """Fingerprint available OCR models, without OCR or remote lookups.

    With no language argument this identifies all installed models: the actual
    per-document language is resolved later by scan detection and BibTeX pins.
    """
    result = {"policy": NATIVE_TEXT_RECOVERY_POLICY, "dpi": list(_DPIS), "psm": 6,
              "max_regions_per_page": _MAX_REGIONS_PER_PAGE, "timeout_seconds": _TIMEOUT,
              "max_regions_per_document": _MAX_REGIONS_PER_DOCUMENT,
              "available": False, "models": {},
              "pymupdf_version": importlib.metadata.version('pymupdf')}
    executable = shutil.which('tesseract')
    if not executable:
        return result
    version,directory = _installation(executable,Path(executable).stat().st_mtime_ns,
                                     os.environ.get('TESSDATA_PREFIX'))
    result['version'] = version
    if not directory or not directory.is_dir():
        return result
    names = sorted(set(languages)) if languages is not None else sorted(p.stem for p in directory.glob('*.traineddata'))
    for name in names:
        path = directory / (name+'.traineddata')
        try:
            stat = path.stat()
            result['models'][name] = _model_hash(str(path),stat.st_mtime_ns,stat.st_size)
        except OSError:
            result['models'][name] = None
    result['available'] = bool(names) and all(result['models'].values())
    return result


def decoded_hint(token):
    """A candidate only; a literal Ã or a printed byte example is not decoded."""
    if not _DAMAGE.search(token):
        return None
    try:
        decoded = token.encode('latin-1').decode('utf-8')
    except UnicodeError:
        return None
    return decoded if (decoded != token and sum(c.isalpha() for c in decoded) >= 2
                       and any(unicodedata.combining(c) for c in unicodedata.normalize('NFD',decoded))) else None


def _clean_token(value):
    return value.strip('.,;:!?()[]{}\"“”„«»')


def supported_reading(hint, readings):
    """Require both raster readings to agree on the complete accented token."""
    hint = _clean_token(hint)
    choices = []
    for reading in readings:
        matches = set()
        for word in reading.split():
            word = _clean_token(word)
            # A damaged text layer may omit a short ending. Every observed
            # decoded character, accent and digit must survive in order.
            if word.casefold().startswith(hint.casefold()) and 0 <= len(word)-len(hint) <= 3:
                if re.findall(r'\d+',word) == re.findall(r'\d+',hint):
                    matches.add(word)
        choices.append(matches)
    if len(choices) != len(_DPIS):
        return None
    common = set.intersection(*choices)
    return next(iter(common)) if len(common) == 1 else None


def _page_regions(page):
    import fitz
    words = page.get_text('words')
    grouped = defaultdict(list)
    for word in words:
        raw = word[4]
        hint = decoded_hint(raw)
        if hint:
            box = fitz.Rect(word[:4])*page.rotation_matrix
            grouped[(word[5],word[6])].append({'original': raw,'decoded_hint':hint,'bbox':list(box)})
    lines = {}
    for block_no,block in enumerate(page.get_text('dict').get('blocks',[])):
        for line_no,line in enumerate(block.get('lines',[])):
            lines[(block_no,line_no)] = fitz.Rect(line['bbox'])*page.rotation_matrix
    regions = []
    for key,candidates in grouped.items():
        rect = lines.get(key)
        if rect is None:
            rect = fitz.Rect(candidates[0]['bbox'])
            for candidate in candidates[1:]:
                rect |= fitz.Rect(candidate['bbox'])
        # Native/OCR font boxes can lie several pixels above actual ink. A
        # neighboring line supplies enough context; one tight word crop does
        # not (the source pilot's PSM7 reading was empty).
        pad = max(12,rect.height*1.2)
        rect = (rect+(-3,-pad,3,pad)) & page.rect
        regions.append({'bbox':list(rect),'candidates':candidates})
    return regions


def inspect_native_text_regions(pdf_path, languages):
    """Read original source candidates before full-page OCR can erase them."""
    import fitz
    languages = list(languages or [])
    report = {'method':NATIVE_TEXT_RECOVERY_POLICY,'languages':languages,'regions':[],
              'candidate_count':0,'confirmed_count':0,'unresolved':[]}
    try:
        pdf = fitz.open(pdf_path)
    except (OSError,fitz.FileDataError):
        report['status'] = 'source_unavailable'
        return report
    executable = shutil.which('tesseract')
    inspected_regions = 0
    with pdf:
        for page_no,page in enumerate(pdf,1):
            regions = _page_regions(page)
            report['candidate_count'] += sum(len(r['candidates']) for r in regions)
            for index,region in enumerate(regions):
                entry = {'page':page_no,'page_size':[page.rect.width,page.rect.height],**region}
                budget_exceeded = index >= _MAX_REGIONS_PER_PAGE or inspected_regions >= _MAX_REGIONS_PER_DOCUMENT
                if budget_exceeded or not executable or not languages:
                    reason = 'region_budget_exceeded' if budget_exceeded else 'ocr_unavailable'
                    report['unresolved'].append({'page':page_no,'bbox':region['bbox'],'reason':reason,
                                                 'original_tokens':[c['original'] for c in region['candidates']]})
                    continue
                inspected_regions += 1
                readings = []
                for dpi in _DPIS:
                    pix = page.get_pixmap(clip=fitz.Rect(region['bbox']),dpi=dpi)
                    try:
                        result = subprocess.run([executable,'stdin','stdout','-l','+'.join(languages),'--psm','6'],
                                                input=pix.tobytes('png'),capture_output=True,timeout=_TIMEOUT)
                        readings.append(result.stdout.decode('utf-8',errors='replace') if result.returncode == 0 else '')
                    except (OSError,subprocess.TimeoutExpired):
                        readings.append('')
                for candidate in entry['candidates']:
                    corrected = supported_reading(candidate['decoded_hint'],readings)
                    candidate['confirmed_text'] = corrected
                    candidate['status'] = 'confirmed' if corrected else 'unresolved'
                    if corrected:
                        report['confirmed_count'] += 1
                    else:
                        report['unresolved'].append({'page':page_no,'bbox':candidate['bbox'],
                            'original':candidate['original'],'reason':'regional_ocr_did_not_agree_with_decoded_accent'})
                report['regions'].append(entry)
    if report['candidate_count']:
        report['producer'] = native_text_recovery_producer(languages)
        with Path(pdf_path).open('rb') as stream:
            report['source_pdf_sha256'] = hashlib.file_digest(stream,'sha256').hexdigest()
    return report


def _nearby_word(page, candidate):
    import fitz
    original_box = fitz.Rect(candidate['bbox'])
    choices = []
    for word in page.get_text('words'):
        box = fitz.Rect(word[:4])*page.rotation_matrix
        intersection = original_box & box
        if intersection.is_empty or intersection.get_area() < .45*min(original_box.get_area(),box.get_area()):
            continue
        value = _clean_token(word[4])
        expected = _clean_token(candidate['confirmed_text'])
        if re.findall(r'\d+',value) != re.findall(r'\d+',expected):
            continue
        if SequenceMatcher(None,value.casefold(),expected.casefold()).ratio() >= .72:
            choices.append((value,list(box)))
    return choices[0] if len(choices) == 1 else None


def apply_native_text_recovery(document, pdf_path, recovery):
    """Apply confirmed original-region words to their prepared-text locations."""
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta, FloatingMeta
    from .source_layout import item_bounds
    report = {'method':NATIVE_TEXT_RECOVERY_POLICY,'repairs':[],'unresolved':[]}
    if not recovery or recovery.get('method') != NATIVE_TEXT_RECOVERY_POLICY:
        return report
    report['producer'] = recovery.get('producer')
    report['source_pdf_sha256'] = recovery.get('source_pdf_sha256')
    report['unresolved'].extend(recovery.get('unresolved',[]))
    with fitz.open(pdf_path) as pdf:
        for region in recovery.get('regions',[]):
            page_no = region['page']
            if not 1 <= page_no <= len(pdf):
                continue
            page = pdf[page_no-1]
            if any(abs(a-b) > .5 for a,b in zip(region['page_size'],[page.rect.width,page.rect.height])):
                report['unresolved'].append({'page':page_no,'reason':'prepared_page_geometry_changed'})
                continue
            for candidate in region['candidates']:
                if not candidate.get('confirmed_text'):
                    continue
                located = _nearby_word(page,candidate)
                if located is None:
                    report['unresolved'].append({'page':page_no,'original':candidate['original'],
                                                 'bbox':candidate['bbox'],
                                                 'reason':'prepared_word_geometry_or_text_unconfirmed'})
                    continue
                old,word_box = located
                corrected = _clean_token(candidate['confirmed_text'])
                owners = []
                for item in document.texts:
                    box = item_bounds(item,document)
                    if box and box[0] == page_no:
                        intersection = fitz.Rect(box[1:]) & fitz.Rect(word_box)
                        if not intersection.is_empty and intersection.get_area() >= .5*fitz.Rect(word_box).get_area():
                            owners.append((item,item))
                for table in document.tables:
                    if len(table.prov) == 1 and table.prov[0].page_no == page_no:
                        for cell in table.data.table_cells:
                            if cell.bbox:
                                b=cell.bbox.to_top_left_origin(page.rect.height)
                                intersection=fitz.Rect(b.l,b.t,b.r,b.b)&fitz.Rect(word_box)
                                if not intersection.is_empty and intersection.get_area() >= .5*fitz.Rect(word_box).get_area():
                                    owners.append((cell,table))
                matching = []
                for item,owner in owners:
                    pattern = re.compile(r'(?<!\w)'+re.escape(old)+r'(?!\w)')
                    matches=list(pattern.finditer(item.text))
                    if len(matches)==1:
                        matching.append((item,owner,matches[0]))
                    elif corrected in item.text:
                        pass  # Already applied to this structured item.
                if len(matching)!=1:
                    if not any(corrected in item.text for item,_ in owners):
                        report['unresolved'].append({'page':page_no,'original':candidate['original'],
                                                     'bbox':candidate['bbox'],
                                                     'reason':'structured_token_alignment_unconfirmed'})
                    continue
                item,owner,match=matching[0]
                if match[0] == corrected:
                    continue
                note={'item_ref':owner.self_ref,'page':page_no,'charspan':list(match.span()),
                      'original':match[0],'replacement':corrected,'original_native_token':candidate['original'],
                      'decoded_hint':candidate['decoded_hint'],'source_bbox':candidate['bbox'],
                      'prepared_bbox':word_box,'evidence':'two_regional_raster_readings_and_prepared_geometry',
                      'dpi':list(_DPIS),'languages':recovery.get('languages',[]),
                      'source_pdf_sha256':recovery.get('source_pdf_sha256'),
                      'method':NATIVE_TEXT_RECOVERY_POLICY,'status':'repaired'}
                item.text=item.text[:match.start()]+corrected+item.text[match.end():]
                report['repairs'].append(note)
                if owner.meta is None:
                    owner.meta=BaseMeta() if owner is item else FloatingMeta()
                notes=list(getattr(owner.meta,'corpus__native_text_recovery',[]) or [])
                if note not in notes:
                    notes.append(note)
                owner.meta.corpus__native_text_recovery=notes
        for problem in report['unresolved']:
            source_box = problem.get('bbox')
            if not source_box:
                continue
            for item in [*document.texts,*document.tables]:
                box = item_bounds(item,document)
                if not box or box[0] != problem.get('page'):
                    continue
                intersection = fitz.Rect(box[1:]) & fitz.Rect(source_box)
                if intersection.is_empty or intersection.get_area() < .5*fitz.Rect(source_box).get_area():
                    continue
                note = {'item_ref':item.self_ref,'status':'unresolved',
                        'method':NATIVE_TEXT_RECOVERY_POLICY,**problem}
                if item.meta is None:
                    item.meta = FloatingMeta() if item in document.tables else BaseMeta()
                notes=list(getattr(item.meta,'corpus__native_text_recovery',[]) or [])
                if note not in notes:
                    notes.append(note)
                item.meta.corpus__native_text_recovery=notes
    return report

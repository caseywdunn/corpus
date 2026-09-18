"""Frozen source-transcription-only proposal population; never reads build outputs."""
import collections,datetime,hashlib,json,re
from pathlib import Path
ROOT=Path('/home/claude/repos/siphonophores')
OUT=Path('/tmp/corpus-v15-notation-holdout')
SEED='20260918-issue303-source-only-v1'
EXCLUDE={'Sutherland_etal2019b','Haddock_etal2005','Hauss_et2016','Pakhomov_etal2000','Kidwai_Amjad2000',
         'Ahuja_etal2026','Bernstein1934','Hosiaetal2024'}
PATTERNS={
 'uncertainty':r'\d[\d.,]*\s*[±∓]\s*\d[\d.,]*',
 'exponent_or_inverse_unit':r'(?:\b(?:[cmknsug]+|sec|min|kg|mol|10)\s*(?:\^\s*[-−]?\d+|[⁻⁺][⁰¹²³⁴⁵⁶⁷⁸⁹]+|[−-][1-9])(?!\d)|\b(?:m|cm|mm)[²³])',
 'arithmetic_or_relation':r'(?:\d[\d.,]*\s*[+×=<>≤≥−]\s*\d[\d.,]*|[A-Za-zα-ω]\s*[=<>≤≥]\s*[-−]?\d[\d.,]*)',
 'measurement_unit':r'\d[\d.,]*\s*(?:[µμ]m|mm|cm|nm|°\s*C|sec|m)(?![\w²³])',
 'ordinary_range':r'(?<!\w)\d[\d.,]*\s*[–—-]\s*\d[\d.,]*(?!\w)',
 'superscript_or_footnote_control':r'(?:[A-Za-zÀ-ö][¹²³⁴⁵⁶⁷⁸⁹⁰*†‡]+|(?<!\w)[¹²³⁴⁵⁶⁷⁸⁹⁰]+(?!\w)|\[FOOTNOTE\]\s*[*†‡¹²³⁴⁵⁶⁷⁸⁹⁰0-9]+)',
}
def digest(data):return hashlib.sha256(data).hexdigest()
def main():
 manifest_path=ROOT/'transcriptions/sources.json'
 manifest=json.loads(manifest_path.read_text())
 policy={'seed':SEED,'source_manifest':str(manifest_path),'source_manifest_sha256':digest(manifest_path.read_bytes()),
 'allowed_inputs':['sources.json','TRANSCRIPTION_PROTOCOL.md','CROSSCHECK_REPORT.md','eligible page_NNN.txt independent transcriptions','original source PDF raster only'],
 'forbidden_inputs':['candidate build/extraction/chunk/embedding/repair outputs','native PDF text or OCR during label verification'],
 'excluded_papers':sorted(EXCLUDE),'exclusion_reason':'Named #303 regressions and gold papers whose candidate outputs were already discussed or inspected in this thread.',
 'max_expressions':24,'initial_stratum_quota':4,'max_per_document':3,'max_per_page':1,'patterns':PATTERNS,
 'ranking':'ascending SHA256(seed|stratum|stem|physical_page|transcription_char_start); one pass per stratum in listed order; then remaining candidates by rank subject to caps until24',
 'deduplication':'Each pattern match is an occurrence. Discard matches whose start falls within40 characters of an earlier-priority pattern match on same page. Select at most one occurrence per page.',
 'verification':'Render selected physical pages from hashed originals. Classify actual scientific meaning or negative; retain all selected unknowns/rejections without swapping after visual review. Transcription guidance is not final source truth.',
 'scope':'Broader source-selected bounded evaluation population; intentionally stratified and capped, not representative of the corpus; final recall/precision require later candidate comparison.',
 'predeclared_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 policy_path=OUT/'selection-policy.json'
 if policy_path.exists():raise RuntimeError('refusing to overwrite frozen policy')
 policy_path.write_text(json.dumps(policy,ensure_ascii=False,indent=2)+'\n')
 # Freeze the policy before opening any candidate transcription text.
 rows=[];transcription_hashes={}
 for stem,source in manifest.items():
  if stem in EXCLUDE:continue
  for page_file in sorted((ROOT/'transcriptions'/stem).glob('page_*.txt')):
   data=page_file.read_bytes();text=data.decode('utf-8');page=int(page_file.stem.split('_')[-1])
   transcription_hashes[str(page_file.relative_to(ROOT))]=digest(data)
   used=[]
   for kind,pattern in PATTERNS.items():
    for match in re.finditer(pattern,text):
     if any(abs(match.start()-pos)<40 for pos in used):continue
     used.append(match.start())
     key=f'{SEED}|{kind}|{stem}|{page}|{match.start()}'
     rows.append({'id':digest(key.encode())[:12],'stratum':kind,'stem':stem,'physical_page':page,'transcription':str(page_file.relative_to(ROOT)),
      'transcription_sha256':digest(data),'charspan':[match.start(),match.end()],'transcribed_expression':match.group(),
      'context':text[max(0,match.start()-120):min(len(text),match.end()+160)],'source_pdf':source['pdf'],'source_sha256':source['sha256'],'rank':digest(key.encode())})
 selected=[];papers=collections.Counter();pages=set()
 def select(row):
  if len(selected)>=24 or papers[row['stem']]>=3 or (row['stem'],row['physical_page']) in pages:return False
  selected.append(row);papers[row['stem']]+=1;pages.add((row['stem'],row['physical_page']));return True
 for kind in PATTERNS:
  picked=0
  for row in sorted((r for r in rows if r['stratum']==kind),key=lambda r:r['rank']):
   if picked==4:break
   if select(row):picked+=1
 for row in sorted(rows,key=lambda r:r['rank']):select(row)
 for index,row in enumerate(selected,1):row['sample_number']=index
 (OUT/'transcription-inputs.json').write_text(json.dumps(transcription_hashes,indent=2)+'\n')
 (OUT/'population.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
 (OUT/'selected.json').write_text(json.dumps(selected,ensure_ascii=False,indent=2)+'\n')
 receipt={'policy_sha256':digest(policy_path.read_bytes()),'selection_script_sha256':digest(Path(__file__).read_bytes()),
  'eligible_transcription_pages':len(transcription_hashes),'candidate_occurrences':len(rows),'population_by_stratum':dict(collections.Counter(r['stratum']for r in rows)),
  'selected_count':len(selected),'selected_by_stratum':dict(collections.Counter(r['stratum']for r in selected)),
  'selected_by_document':dict(papers),'selected_sha256':digest((OUT/'selected.json').read_bytes()),
  'frozen_before_visual_review_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 (OUT/'selection-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))
 for row in selected:print(row['sample_number'],row['stem'],row['physical_page'],row['stratum'],repr(row['transcribed_expression']),repr(row['context']))
if __name__=='__main__':main()

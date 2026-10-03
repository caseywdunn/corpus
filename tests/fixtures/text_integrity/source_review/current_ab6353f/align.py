import json,re,sys
from pathlib import Path
P=Path(sys.argv[1]);labels=json.loads(Path('/home/claude/repos/corpus/tests/fixtures/text_integrity/source_review/labels.json').read_text());states={r['hash']:r for r in json.loads((P/'initial-state-receipt.json').read_text())['documents']}
output=[]
for row in labels:
 h=row['source_sha256'][:12];stage=(states[h].get('pipeline_state')or{}).get('stages',{});n=row['sample_number'];entry={'sample_number':n,'stem':row['stem'],'physical_page':row['physical_page'],'selected_source_literal':row['review']['source_literal']}
 if not row['review']['eligible_for_literal_expression_comparison']:
  entry['status']='unscorable_source';output.append(entry);continue
 if not stage.get('text_chunking',{}).get('completed_at'):
  entry['status']='pending';output.append(entry);continue
 d=P/'captured'/h;doc=json.loads((d/'docling_doc.json').read_text());chunks=json.loads((d/'chunks.json').read_text())['chunks'];text=json.loads((d/'text.json').read_text())
 keep=stage['text_chunking']['input_fingerprint'].get('keeppages');entry['keeppages']=keep
 if keep:
  selected_pages=[]
  for part in keep.split(','):
   ends=re.split(r'\s*--?\s*',part.strip());selected_pages.extend(range(int(ends[0]),int(ends[-1])+1))
  if row['physical_page']not in selected_pages:
   entry['status']='excluded_by_source_page_selection';output.append(entry);continue
  page=selected_pages.index(row['physical_page'])+1
 else:page=row['physical_page']
 entry['prepared_page']=page;entry['status']='inspection';size=doc['pages'][str(page)]['size'];W,H=size['width'],size['height'];box=row['review']['rendered_source']['relative_page_bbox'];source_box=[box[0]*W,box[1]*H,box[2]*W,box[3]*H];entry['candidate_size']=[W,H];entry['candidate_region']=source_box
 selected=[]
 for category in ['texts','tables','pictures']:
  for item in doc[category]:
   for prov in item.get('prov',[]):
    if prov['page_no']!=page:continue
    b=prov['bbox'];bounds=[b['l'],H-b['t'],b['r'],H-b['b']]if b['coord_origin']=='BOTTOMLEFT'else[b['l'],b['t'],b['r'],b['b']]
    intersects=min(bounds[2],source_box[2])>max(bounds[0],source_box[0])and min(bounds[3],source_box[3])>max(bounds[1],source_box[1])
    if intersects:
     selected.append({'category':category,'item_ref':item['self_ref'],'label':item['label'],'text':item.get('text'),'orig':item.get('orig'),'meta':item.get('meta'),'bbox_top_left':bounds,'provenance':prov,'data':item.get('data')if category=='tables'else None})
 refs={s['item_ref']for s in selected};entry['region_items']=selected;entry['linked_chunks']=[c for c in chunks if any(s.get('item_ref')in refs for s in c.get('source_items',[]))];report=text['source_text_integrity']['scientific_notation'];entry['page_scientific_repairs']=[r for r in report['repairs']if r.get('page')==page];entry['page_scientific_unresolved']=[r for r in report['unresolved']if r.get('page')==page];entry['total_scientific_repairs_in_document']=len(report['repairs']);output.append(entry)
(P/'region-candidates.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
for e in output:
 print('\nSAMPLE',e['sample_number'],e['stem'],e['status'],'physical',e['physical_page'],'prepared',e.get('prepared_page'),'keep',e.get('keeppages'))
 for item in e.get('region_items',[]):
  print(item['item_ref'],item['label'],str(item.get('text')or item.get('data'))[:1800])
 print('linked chunks',[c['chunk_id']for c in e.get('linked_chunks',[])],'page repair count',len(e.get('page_scientific_repairs',[])))

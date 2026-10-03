"""Read-only current-gold capture and conservative replay of frozen source grades."""
import argparse
import collections
import copy
import datetime
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent
SOURCE=Path('/home/claude/repos/corpus/tests/fixtures/text_integrity/source_review')
BASE=Path('/tmp/corpus-v15-notation-comparison-complete-20260922')
def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()
def digest_obj(value):
    return digest_bytes(json.dumps(value,sort_keys=True).encode())
def load(path):
    return json.loads(path.read_text())
def save(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def complete_current(state,requirements):
    stages=state.get('stages',{})
    for name,expected in requirements['stage_config_sha256'].items():
        stage=stages.get(name,{})
        if not stage.get('completed_at'):
            return False,f'{name}_not_completed'
        if digest_obj(stage.get('input_fingerprint',{}).get('config',{}))!=expected:
            return False,f'{name}_prior_or_different_producer'
    return True,'completed_current_producers'
def projection(region,markdown):
    # Exact evidence equality only. Altered text/layout/chunk boundaries require
    # a new source judgment; current source receipts never inherit old repairs.
    return {'region':region.get('candidate_region'),
            'items':[{k:i.get(k)for k in ['category','item_ref','label','text','bbox_top_left','provenance','data']}for i in region.get('region_items',[])],
            'chunks':[{k:c.get(k)for k in ['chunk_id','text','source_items']}for c in region.get('linked_chunks',[])],
            'markdown_sha256':digest_bytes(markdown.encode())}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gold',type=Path,default=Path('/tmp/corpus-v15-acceptance/gold'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--require-complete',action='store_true',help='Refuse final mode unless all35 gold documents have matching completed current stages.')
    parser.add_argument('--adjudications',type=Path,default=SCRIPT/'adjudications.json')
    args=parser.parse_args();out=args.out.resolve()
    if out.exists():parser.error('Output must be a new scratch directory; historical snapshots are never overwritten.')
    if out==args.gold.resolve()or args.gold.resolve()in out.parents:parser.error('Output must be outside the gold input/output tree.')
    requirements=load(SCRIPT/'producer-requirements.json')
    for filename,key in [('refresh-20260922-input.json','refresh_input_sha256'),('config-v15-sep22.yaml','config_sha256'),('input-receipt.json','gold_input_receipt_sha256')]:
        assert digest_bytes((args.gold/filename).read_bytes())==requirements[key],filename
    assert digest_bytes((SOURCE/'labels.json').read_bytes())==requirements['labels_sha256']
    labels=load(SOURCE/'labels.json');input_receipt=load(args.gold/'input-receipt.json')
    inputs=[r for r in input_receipt['inputs']if r.get('path','').lower().endswith('.pdf')]
    assert len(inputs)==input_receipt['gold_document_count']
    indexed={r['sha256']:r for r in inputs}
    assert all(r['source_sha256']in indexed for r in labels)
    population=[]
    for item in inputs:
        d=args.gold/'output/documents'/item['sha256'][:12]
        try:state=load(d/'pipeline_state.json')
        except (FileNotFoundError,json.JSONDecodeError):state={}
        ready,reason=complete_current(state,requirements)
        population.append({'stem':item['stem'],'source_sha256':item['sha256'],'completed_current':ready,'reason':reason})
    out.mkdir(parents=True)
    save(out/'population-state.json',population)
    now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    save(out/'run-receipt.json',{'captured_utc':now,'code_revision':requirements['code_revision'],'scope':'completed current-gold sample'if args.require_complete else'provisional completed-current-document snapshot','producer_requirements':requirements,'script_sha256':digest_bytes(Path(__file__).read_bytes()),'gold_population_complete':all(r['completed_current']for r in population),'source_labels_unchanged':True,'no_new_pipeline_or_models':True})
    if args.require_complete and not all(r['completed_current']for r in population):
        print('FINAL REFUSED:',sum(r['completed_current']for r in population),'of',len(population),'gold documents have matching completed current stages; state receipt saved.',file=sys.stderr)
        return 2
    states=[];captures=[]
    for row in {r['source_sha256']:r for r in labels}.values():
        h=row['source_sha256'][:12];d=args.gold/'output/documents'/h
        try:before=(d/'pipeline_state.json').read_bytes();state=json.loads(before)
        except (FileNotFoundError,json.JSONDecodeError):before=b'';state={}
        ready,reason=complete_current(state,requirements)
        record={'stem':row['stem'],'hash':h,'source_sha256':row['source_sha256'],'completed_current':ready,'reason':reason}
        if ready:
            dest=out/'captured'/h;dest.mkdir(parents=True);files={}
            try:
                for name in ['pipeline_state.json','docling_doc.json','text.json','chunks.json','figures.json','metadata.json','scan_detection.json','summary.json']:
                    src=d/name
                    if src.exists():
                        data=src.read_bytes();(dest/name).write_bytes(data);files[name]={'sha256':digest_bytes(data),'size':len(data)}
                assert before==(d/'pipeline_state.json').read_bytes()
                assert all(n in files for n in ['docling_doc.json','text.json','chunks.json'])
            except (AssertionError,FileNotFoundError,json.JSONDecodeError):
                shutil.rmtree(dest);record.update(completed_current=False,reason='changed_or_incomplete_during_capture');state={}
            else:record.update(files=files,state_stable_during_read=True);captures.append(record)
        else:state={}  # Old completed stages cannot pass the inherited aligner.
        states.append({**record,'pipeline_state':state})
    save(out/'initial-state-receipt.json',{'captured_utc':now,'documents':states})
    save(out/'capture-receipt.json',captures)
    with (out/'alignment.txt').open('w')as log:
        subprocess.run([sys.executable,str(SCRIPT/'align.py'),str(out)],stdout=log,check=True)
    regions=load(out/'region-candidates.json');baseline_regions={r['sample_number']:r for r in load(BASE/'region-candidates.json')};baseline_grades={r['sample_number']:r for r in load(BASE/'results.json')}
    captures_by_hash={r['hash']:r for r in captures};adjudications=load(args.adjudications)if args.adjudications.exists()else{}
    results=[]
    for label,region in zip(labels,regions,strict=True):
        n=label['sample_number'];h=label['source_sha256'][:12]
        row={'sample_number':n,'stem':label['stem'],'source_sha256':label['source_sha256'],'physical_page':label['physical_page'],'prepared_page':region.get('prepared_page'),'source_literal':label['review']['source_literal'],'status':region['status']}
        if region['status']in ['pending','excluded_by_source_page_selection','unscorable_source']:
            row['surface_status']={k:region['status']for k in ['docling','markdown','chunks']}
            if region.get('keeppages'):row['keeppages']=region['keeppages']
        else:
            files=captures_by_hash[h]['files'];key=digest_obj({k:files[k]['sha256']for k in ['docling_doc.json','text.json','chunks.json','figures.json']})
            row['candidate_artifact_key']=key
            current_markdown=load(out/'captured'/h/'text.json')['text'];old_markdown=load(BASE/'captured'/h/'text.json')['text']
            if projection(region,current_markdown)==projection(baseline_regions[n],old_markdown):
                grade=baseline_grades[n]
                row.update(status='graded_by_exact_unchanged_surface_evidence',surface_status={k:grade[f'{v}_status']for k,v in [('docling','docling'),('markdown','saved_markdown'),('chunks','chunks')]},assessment=grade['assessment'],candidate_item_refs=grade.get('candidate_item_refs'),candidate_chunk_ids=grade.get('candidate_chunk_ids'))
            elif str(n)in adjudications and adjudications[str(n)]['candidate_artifact_key']==key:
                row.update(copy.deepcopy(adjudications[str(n)]));row['status']='source_reviewed_current'
            else:
                row.update(status='review_required_changed_surface_evidence',surface_status={k:'review_required'for k in ['docling','markdown','chunks']})
            if n==22:
                pattern=r'750\s*[-–—]\s*500\s*m\.?'
                region_hits=[{'item_ref':i['item_ref'],'bbox':i['bbox_top_left'],'text':i.get('text')}for i in region['region_items']if re.search(pattern,i.get('text')or'')]
                chunk_hits=[{'chunk_id':c['chunk_id'],'text':c['text']}for c in region['linked_chunks']if re.search(pattern,c['text'])]
                row['selected_table_cell_probe']={'source_region_item_hits':region_hits,'linked_chunk_hits':chunk_hits,'markdown_hits':[{'charspan':list(m.span()),'context':current_markdown[max(0,m.start()-100):m.end()+100]}for m in re.finditer(pattern,current_markdown)],'warning':'Probe is source-region/context evidence, not automatic table-structure or semantic grading. Review neighboring station1772/netN7OV before credit.'}
        results.append(row)
    save(out/'results.json',results)
    summary={'scope':'Current producer snapshot only; no baseline rows counted as current','sample_size':24,'counts':{surface:dict(collections.Counter(r['surface_status'][surface]for r in results))for surface in ['docling','markdown','chunks']},'population_completed_current':sum(r['completed_current']for r in population),'population_total':len(population),'current_sample_documents_captured':len(captures),'review_required_samples':[r['sample_number']for r in results if r['status']=='review_required_changed_surface_evidence'],'pending_samples':[r['sample_number']for r in results if r['status']=='pending'],'repair_precision':'Not inherited from prior-build operations; inspect current source-aligned receipts separately.','final_fully_adjudicated':not any(r['status']in ['pending','review_required_changed_surface_evidence']for r in results)}
    save(out/'summary.json',summary);print(json.dumps(summary,indent=2))
    if args.require_complete and not summary['final_fully_adjudicated']:
        print('Artifacts are complete but changed regions require explicit source adjudication; this is not a final graded score.',file=sys.stderr)
        return 3
    return 0
if __name__=='__main__':raise SystemExit(main())

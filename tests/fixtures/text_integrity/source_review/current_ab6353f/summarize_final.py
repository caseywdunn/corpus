"""Summarize already source-adjudicated captures; never grades source images."""
import argparse, collections, hashlib, json
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('capture',type=Path)
a=p.parse_args(); root=a.capture
load=lambda path:json.loads(path.read_text())
save=lambda path,value:path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
rows=load(root/'results.json'); basic=load(root/'summary.json'); run=load(root/'run-receipt.json')
assert run['gold_population_complete'] and basic['final_fully_adjudicated']
assert len(rows)==24 and all(r['status'] not in ('pending','review_required_changed_surface_evidence') for r in rows)
by_num={r['sample_number']:r for r in rows}
regions={r['sample_number']:r for r in load(root/'region-candidates.json')}
cap={r['hash']:r for r in load(root/'capture-receipt.json')}
counts={s:dict(collections.Counter(r['surface_status'][s] for r in rows)) for s in ['docling','markdown','chunks']}
for c in counts.values():
 assert c.get('excluded_by_source_page_selection')==4 and c.get('unscorable_source')==1
 assert sum(c.values())==24
assert counts['markdown']==counts['chunks']

operations=[]
for n,item,span in [(1,'#/texts/126',[30,33]),(2,'#/texts/67',[347,350])]:
 row=by_num[n];h=row['source_sha256'][:12]
 report=load(root/'captured'/h/'text.json')['source_text_integrity']['scientific_notation']
 hits=[r for r in report['repairs'] if r['item_ref']==item and r.get('charspan')==span]
 assert len(hits)==1 and hits[0]['replacement']=='⁻³'
 assert row['surface_status']==dict.fromkeys(['docling','markdown','chunks'],'faithful')
 operations.append({'sample_number':n,'source_sha256':row['source_sha256'],'physical_page':row['physical_page'],'method':report['method'],'receipt':hits[0],'text_artifact_sha256':cap[h]['files']['text.json']['sha256'],'source_review':'Current selected operation matches the frozen original source crop; independent selected-operation review supplied separately.','source_correct':True})
row=by_num[4];h=row['source_sha256'][:12];report=load(root/'captured'/h/'text.json')['source_text_integrity']['pdf_cmap']
macro=[r for r in report['repairs'] if r['item_ref']=='#/texts/28'];assert len(macro)==1
entry=[r for r in macro[0]['formatting'] if r['char_offset']==1249];assert len(entry)==1 and entry[0]['replacement']=='³'
assert all(x['text'].strip()=='3' for x in entry[0]['raster']['observations'])
operations.append({'sample_number':4,'source_sha256':row['source_sha256'],'physical_page':row['physical_page'],'method':report['method'],'producer':report['producer'],'item_ref':'#/texts/28','receipt':entry[0],'text_artifact_sha256':cap[h]['files']['text.json']['sha256'],'scope':'Only the first/surface mg/m³ raised-digit formatting suboperation. Whole-paragraph CMap mappings and the middle/bottom unit formatting entries are excluded.','source_review':'Frozen source crop04 and exact current formatting position independently reviewed.','source_correct':True})
save(root/'selected-current-operations.json',{'operations':operations,'strata':{'source_glyph_alignment_v1':{'correct':2,'reviewed':2,'papers':1},'exact_regional_type1_scalar_cmap_raised_digit':{'correct':1,'reviewed':1,'papers':1}},'aggregation':'Three source-correct selected atomic exponent operations across two different producers and two papers. Not complete paragraph/CMap precision, proposal precision, or corpus-wide precision. No historical operation counts inherited.'})

# Current refusal is preserved independently of the neighboring repaired2000.
r=by_num[16];h=r['source_sha256'][:12];t=load(root/'captured'/h/'text.json')
original=t['source_text_integrity']['original_native_layer']
def walk(v):
 if isinstance(v,dict):
  if v.get('reason')=='digit_ocr_disagreement' and v.get('source_evidence',{}).get('prefix')=='0.3-10мм':yield v
  for x in v.values():yield from walk(x)
 elif isinstance(v,list):
  for x in v:yield from walk(x)
refusals=list(walk(original));assert len(refusals)==1
save(root/'selected-current-unresolved.json',{'sample_number':16,'surface_status':r['surface_status'],'item_ref':'#/texts/125','chunk_id':'chunk_97','text_artifact_sha256':cap[h]['files']['text.json']['sha256'],'selected_receipt':refusals[0],'not_credited':'Neighboring2000мм³ belongs to another occurrence.'})

metrics={}
for surf,c in counts.items():
 f=c.get('faithful',0);bad=c.get('corrupted',0);unknown=c.get('fidelity_indeterminate',0);emitted=f+bad+unknown
 metrics[surf]={'expression_fidelity_adjudicated':f'{f}/{f+bad}','expression_fidelity_all_emitted_bounds':[f'{f}/{emitted}',f'{f+unknown}/{emitted}'],'source_expression_recall_consumed_page_bounds':[f'{f}/19',f'{f+unknown}/19'],'source_expression_recall_all_verified_source_bounds':[f'{f}/23',f'{f+unknown}/23']}
 # No targeted notation recovery at the selected expression; OCR/other repairs elsewhere can exist.
 noop=[r for r in rows if r['sample_number']not in [1,2,4,8] and r['surface_status'][surf]in ['faithful','corrupted','fidelity_indeterminate']]
 nc=collections.Counter(r['surface_status'][surf] for r in noop);nf=nc['faithful'];nu=nc['fidelity_indeterminate']
 metrics[surf]['no_selected_notation_repair_emitted_expression_correctness_bounds']=[f'{nf}/{len(noop)}',f'{nf+nu}/{len(noop)}']
strata={}
for name,numbers in [('fix_informed_regression',[4,8,16]),('other_frozen_source_selections',[n for n in range(1,25)if n not in [4,8,16]])]:
 strata[name]={'sample_numbers':numbers,'counts':{s:dict(collections.Counter(by_num[n]['surface_status'][s] for n in numbers))for s in counts}}
summary={'candidate_code_revision':run['code_revision'],'candidate_scope':'Completed normal current CPU extraction/chunking refresh; current stages checked for all35 gold PDFs, all14 sampled documents captured after matching completion. No model or OCR run by the scorer.','source_label_sha256':run['producer_requirements']['labels_sha256'],'denominators':{'selected':24,'source_verified':23,'unscorable_not_printed':1,'excluded_by_recorded_keeppages':4,'verified_on_consumed_pages':19,'pending':0,'review_required':0},'counts':counts,'metrics':metrics,'selected_current_atomic_repair_precision':{'source_glyph_signed_exponent':'2/2','cmap_first_unit_raised_digit':'1/1','combined_descriptive_count':'3/3','papers':2,'scope':'Current source-aligned reviewed atomic operations only; heterogeneous producer strata stay explicit. Whole-paragraph CMap correctness is not inferred.'},'selected_damage_recovery':{'samples':[1,2,4,8,16],'expression_recovery':'4/5','scope':'Five selected expressions with explicit repair-relevant damaged observations. This is not independently enumerated proposal/detector recall.'},'source_strata':strata,'unmeasured':{'expression_output_precision':'No independently exhaustive enumeration of all spurious/contradictory numerical outputs; expression fidelity is not relabeled precision.','proposal_precision':'No separately frozen proposal population, including rejections.','whole_sample_embedding_and_serving':'This receipt stops at saved CPU chunks; a separately authorized verifier owns vectors and live-tool propagation.'},'limits':['Source-first selection is independent of the named regressions, but Chen4/8 and Russian16 subsequently informed fixes. Current results are regression-sample scores, not fresh unseen holdout accuracy.','No selected±/∓,µm,or subtraction equation; all inverse-unit selections come from one Mańko paper.','Four keeppages exclusions remain coverage limits. Figure-axis/footer and any table-to-picture omissions stay in consumed-page recall.','Caption magnification×→x/X and decimal middle-dot→point substitutions are explicit literal discrepancies accepted only in their source-proved caption roles.','Expression correctness does not establish correctness of surrounding Chinese/Japanese prose or whole CMap paragraphs.','Footnote role is not inferred from a baseline digit in markdown/chunks; explicit uncertainty is retained.']}
save(root/'graded-summary.json',summary)
print(json.dumps({'counts':counts,'metrics':metrics,'selected_operations':3,'source_strata':strata},ensure_ascii=False,indent=2))

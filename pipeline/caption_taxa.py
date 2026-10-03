"""Build caption-specific taxon evidence after all figure producers (#323)."""
from __future__ import annotations

import hashlib
import json
import logging
import re
from collections import Counter

from .taxa import TaxonomyDB, _ABBREV_BINOMIAL_RE, _expand_abbreviation

PRODUCER = 'caption-taxa-v1'
logger = logging.getLogger(__name__)
_TOKENS = re.compile(r'[^\W\d_][\w-]*', re.UNICODE)
_ABBREVIATIONS = re.compile(_ABBREV_BINOMIAL_RE.pattern, re.VERBOSE | re.IGNORECASE)


def create_schema(conn):
    conn.executescript('''
        CREATE TABLE IF NOT EXISTS caption_taxon_evidence (
            corpus_hash TEXT NOT NULL, figure_id TEXT NOT NULL, evidence_json TEXT NOT NULL,
            PRIMARY KEY(corpus_hash,figure_id));
        CREATE TABLE IF NOT EXISTS caption_taxon_links (
            corpus_hash TEXT NOT NULL, figure_id TEXT NOT NULL, taxon_id TEXT NOT NULL,
            PRIMARY KEY(corpus_hash,figure_id,taxon_id));
        CREATE INDEX IF NOT EXISTS idx_caption_taxon_links_taxon ON caption_taxon_links(taxon_id);
        CREATE TABLE IF NOT EXISTS caption_taxon_receipts (
            corpus_hash TEXT PRIMARY KEY, input_fingerprint TEXT NOT NULL);
    ''')


def full_names(text, taxonomy, names):
    """Longest exact name spans; case/punctuation never permit partial words."""
    tokens=list(_TOKENS.finditer(text))
    hits=[]
    i=0
    while i<len(tokens):
        found=False
        for size in range(min(4,len(tokens)-i),0,-1):
            start,end=tokens[i].start(),tokens[i+size-1].end()
            surface=text[start:end]
            normalized=' '.join(surface.lower().split())
            if normalized not in names:
                continue
            resolved=taxonomy.lookup(normalized)
            if resolved:
                hits.append({'mention_text':surface,'text_span':[start,end],
                             'matched_name':normalized,'method':'full_name',**resolved})
                if size>1:
                    genus_surface=tokens[i].group()
                    genus=taxonomy.lookup(genus_surface)
                    if genus and (genus.get('rank') or '').lower()=='genus':
                        hits.append({'mention_text':genus_surface,'text_span':[start,tokens[i].end()],
                                     'method':'full_name',**genus})
                i+=size
                found=True
                break
        if not found:
            i+=1
    return hits


def genera(hits):
    return Counter(hit['mention_text'].split()[0].capitalize() for hit in hits)


def caption_evidence(caption, taxonomy, names, document_genera, context_mentions=()):
    matches=full_names(caption,taxonomy,names)
    local=genera(matches)
    unresolved=[]
    for match in _ABBREVIATIONS.finditer(caption):
        start,end=match.span()
        if any(s<end and start<e for s,e in (m['text_span'] for m in matches)):
            continue
        prefix,epithet=match.groups()
        # Caption context is more local than a document's unrelated genera.
        scope='caption' if any(g.lower().startswith(prefix.lower()) for g in local) else 'document'
        context=local if scope=='caption' else document_genera
        expanded,resolved,candidates=_expand_abbreviation(prefix,epithet.lower(),context,taxonomy,names)
        if resolved:
            matches.append({'mention_text':match.group(),'text_span':[start,end],
                'matched_name':expanded,'method':'contextual_abbreviation',
                'context_scope':scope,'context_genera':sorted(context),
                'candidate_names':candidates,**resolved})
            if scope=='document':
                supporting=[m for m in context_mentions
                            if m['mention_text'].split()[0].lower()==expanded.split()[0].lower()]
                matches[-1]['context_mentions']=supporting[:20]
                matches[-1]['context_mentions_total']=len(supporting)
        else:
            reason='ambiguous_abbreviation'
            if not candidates:
                # Diagnostic candidates never count as contextual support.
                possibilities=Counter(name.split()[0].capitalize() for name in names
                                      if ' ' in name and name.split()[1]==epithet.lower())
                _,diagnostic_resolved,candidates=_expand_abbreviation(prefix,epithet.lower(),possibilities,taxonomy,names)
                if diagnostic_resolved:
                    reason='insufficient_genus_context'
            if candidates:
                unresolved.append({'mention_text':match.group(),'text_span':[start,end],
                    'reason':reason,
                    'candidate_names':candidates,'context_scope':scope,'context_genera':sorted(context)})
    matches.sort(key=lambda m:m['text_span'])
    return {'availability':'materialized','matches':matches[:200],
            'matches_total':len(matches),'matches_truncated':len(matches)>200,
            'unresolved':unresolved[:50],'unresolved_total':len(unresolved),
            'unresolved_truncated':len(unresolved)>50,'producer_version':PRODUCER}


def materialize(conn, output_dir):
    """Content receipts cover final figures, document context and taxonomy.

    Parse before replacing a paper's derived rows. Corrupt/missing previously
    indexed figure evidence blocks a successful build, preserving prior rows.
    """
    create_schema(conn)
    docs=output_dir/'documents'
    taxonomy_path=output_dir/'taxonomy.sqlite'
    folders=sorted(p for p in docs.iterdir() if p.is_dir()) if docs.is_dir() else []
    known={sha:fp for sha,fp in conn.execute('SELECT corpus_hash,input_fingerprint FROM caption_taxon_receipts')}
    for sha in known.keys()-{p.name for p in folders}:
        for table in ('caption_taxon_evidence','caption_taxon_links','caption_taxon_receipts'):
            conn.execute(f'DELETE FROM {table} WHERE corpus_hash=?',(sha,))
    stats={'papers':0,'skipped':0,'errors':0}
    if not known and not any((folder/'figures.json').exists() for folder in folders):
        conn.commit()
        return stats
    taxonomy=None
    try:
        taxonomy_sha=hashlib.sha256(taxonomy_path.read_bytes()).hexdigest() if taxonomy_path.exists() else None
        if taxonomy_sha:
            taxonomy=TaxonomyDB(taxonomy_path)
            names=taxonomy.name_set()
        else:
            names=set()
        for folder in folders:
            figure_path=folder/'figures.json'
            if not figure_path.exists():
                if folder.name in known:
                    stats['errors']+=1
                    logger.error('Previously indexed caption evidence is missing: %s',figure_path)
                continue
            try:
                raw=figure_path.read_bytes()
                figures=json.loads(raw)['figures']
                context_path=folder/'chunks.json'
                context_raw=context_path.read_bytes() if context_path.exists() else b'{"chunks":[]}'
                chunks=json.loads(context_raw)['chunks']
                if not isinstance(figures,list) or not isinstance(chunks,list):
                    raise ValueError('Expected figures/chunks lists')
                inputs={'figures_sha256':hashlib.sha256(raw).hexdigest(),
                        'chunks_sha256':hashlib.sha256(context_raw).hexdigest(),
                        'taxonomy_sha256':taxonomy_sha,'producer_version':PRODUCER}
                fingerprint=json.dumps(inputs,sort_keys=True)
                if known.get(folder.name)==fingerprint:
                    stats['skipped']+=1
                    continue
                context=Counter()
                context_mentions=[]
                if taxonomy:
                    for chunk in chunks:
                        chunk_hits=full_names(chunk.get('text') or '',taxonomy,names)
                        context.update(genera(chunk_hits))
                        context_mentions.extend({'chunk_id':chunk.get('chunk_id'),'mention_text':hit['mention_text'],
                                                 'text_span':hit['text_span']} for hit in chunk_hits)
                records=[]
                ids=set()
                for figure in figures:
                    figure_id=figure['figure_id']
                    if not isinstance(figure_id,str) or not figure_id or figure_id in ids:
                        raise ValueError('Missing or duplicate figure ID')
                    ids.add(figure_id)
                    caption=figure.get('caption_text') or figure.get('caption') or ''
                    evidence=caption_evidence(caption,taxonomy,names,context,context_mentions) if taxonomy else {
                        'availability':'no_taxonomy','matches':[],'unresolved':[],'producer_version':PRODUCER}
                    evidence['input_fingerprint']=inputs
                    records.append((figure_id,evidence))
            except (OSError,ValueError,KeyError,TypeError,AttributeError) as exc:
                logger.error('Cannot materialize caption evidence for %s: %s',folder,exc)
                stats['errors']+=1
                continue
            for table in ('caption_taxon_evidence','caption_taxon_links'):
                conn.execute(f'DELETE FROM {table} WHERE corpus_hash=?',(folder.name,))
            for figure_id,evidence in records:
                conn.execute('INSERT INTO caption_taxon_evidence VALUES (?,?,?)',
                             (folder.name,figure_id,json.dumps(evidence,sort_keys=True,ensure_ascii=False)))
                for taxon_id in {str(m['accepted_taxon_id']) for m in evidence['matches']}:
                    conn.execute('INSERT INTO caption_taxon_links VALUES (?,?,?)',(folder.name,figure_id,taxon_id))
            conn.execute('INSERT OR REPLACE INTO caption_taxon_receipts VALUES (?,?)',(folder.name,fingerprint))
            stats['papers']+=1
    finally:
        if taxonomy:
            taxonomy.close()
    conn.commit()
    return stats

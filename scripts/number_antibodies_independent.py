#!/usr/bin/env python3
"""Own offline anchor numbering; public per-position frequency tables only."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
if (ROOT/'src').exists():sys.path.insert(0,str(ROOT/'src'))
from antibody_numbering.engine import AnchorNumberer,SCHEMES


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_fasta(path):
    records=[]
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            if not line[1:].strip():raise ValueError('Empty FASTA identifier')
            records.append({'record_id':f'seq_{len(records)+1:06d}','header':line[1:],'raw_sequence':''})
        else:
            if not records:raise ValueError('Sequence occurs before a FASTA identifier')
            records[-1]['raw_sequence']+=''.join(line.split())
    if not records:raise ValueError('No FASTA input records')
    for r in records:
        raw=r['raw_sequence'].upper()
        r['source_indices']=[i+1 for i,aa in enumerate(raw)if aa!='-']
        r['normalized_sequence']=raw.replace('-','').replace('?','X')
        r['status']=('skipped_stop_character'if'*'in raw else
                     'skipped_invalid_character'if set(raw)-set('ACDEFGHIKLMNPQRSTVWYX?-')else
                     'skipped_empty_sequence'if not r['normalized_sequence']else'ready')
    return records


class FrequencyTables:
    def __init__(self,roots):self.roots=roots;self.cache={};self.used={}
    def lookup(self,scheme,position,aa):
        if not position:return '', 'unmapped_position', ''
        if aa=='X':return '', 'unknown_residue', ''
        key=(scheme,position)
        if key not in self.cache:
            matches=[Path(p)/(scheme+'_'+position+'.json')for p in self.roots]
            matches=[p for p in matches if p.exists()]
            if not matches:self.cache[key]=None
            else:
                # Ordered, explicit precedence; never overwrite a historical table.
                path=matches[0];d=json.loads(path.read_text())
                assert (d['scheme'],d['position'],d['organism'])==(scheme,position,'All')
                assert sum(d['counts'].values())==d['total_including_gap']
                self.cache[key]=(d,path);self.used[str(path)]={'sha256':sha(path),'url':d.get('url'),'fetched_at_utc':d.get('fetched_at_utc')}
        ref=self.cache[key]
        if ref is None:return '', 'missing_frequency_reference', ''
        d,path=ref
        if not d['total_including_gap']:return '', 'empty_reference', str(path)
        return 100*d['counts'][aa]/d['total_including_gap'],'public_position_frequency',str(path)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--frameworks',type=Path,
        default=ROOT/'data/reference/independent_numbering_public_20261109/framework_templates.json',
        help='Public framework JSON; defaults to the bundled project asset: %(default)s')
    p.add_argument('--frequency-cache',action='append',type=Path,default=[])
    args=p.parse_args()
    assert not args.output.exists(),'Use a new output directory'
    records=read_fasta(args.input);numberer=AnchorNumberer(args.frameworks);frequencies=FrequencyTables(args.frequency_cache)
    args.output.mkdir(parents=True)
    columns=['record_id','header','scheme','chain_type','raw_index','normalized_index','input_residue','model_residue','position','status','frequency_percent','frequency_status','frequency_reference','backend']
    summaries=[];predictions={};counts=Counter();missing=set();rowcount=0
    with (args.output/'numbering_with_frequencies.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        for record in records:
            if record['status']!='ready':
                for sc in SCHEMES:summaries.append({'record_id':record['record_id'],'scheme':sc,'status':record['status'],'numbered_residues':0})
                continue
            seq=record['normalized_sequence'];result=numberer.number(seq);predictions[record['record_id']]=result
            for sc in SCHEMES:
                value=result['schemes'][sc]
                if isinstance(value,list):positions=value;status='unclassified'
                else:positions=value['positions'];status=value['status']
                assert len(positions)==len(seq)
                mapped=[x for x in positions if x];assert len(mapped)==len(set(mapped)),(record['record_id'],sc,'duplicate_position')
                counts[sc+':'+status]+=1
                summaries.append({'record_id':record['record_id'],'scheme':sc,'status':status,'numbered_residues':len(mapped)})
                for i,(aa,position)in enumerate(zip(seq,positions)):
                    frequency,fstatus,source=frequencies.lookup(sc,position,aa)
                    if fstatus=='missing_frequency_reference':missing.add((sc,position))
                    writer.writerow({'record_id':record['record_id'],'header':record['header'],'scheme':sc,
                        'chain_type':result.get('chain',''),'raw_index':record['source_indices'][i],'normalized_index':i+1,
                        'input_residue':record['raw_sequence'][record['source_indices'][i]-1],'model_residue':aa,'position':position,
                        'status':status if position or status!='numbered'else'unmapped_position','frequency_percent':frequency,
                        'frequency_status':fstatus,'frequency_reference':source,'backend':'own_anchor_region_v1'})
                    rowcount+=1
    with (args.output/'records_summary.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=['record_id','scheme','status','numbered_residues']);writer.writeheader();writer.writerows(summaries)
    def save(name,value):(args.output/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    save('input_ledger.json',records);save('backend_outputs.json',predictions);save('missing_frequency_positions.json',sorted(missing))
    import inspect
    engine_path=Path(inspect.getfile(AnchorNumberer))
    manifest={'phase':'OFFLINE_NUMBERING_FINISHED_REQUIRES_AUDIT','created_at_utc':datetime.now(timezone.utc).isoformat(),
        'backend':'Own anchor detection, framework alignment and computed five-scheme loop rules',
        'inference_uses_anarci_or_anarcii':False,'inference_uses_query_reference_cache':False,'inference_uses_network':False,
        'public_frameworks':str(args.frameworks),'public_frameworks_sha256':sha(args.frameworks),
        'engine_source_sha256':sha(engine_path),'runner_source_sha256':sha(__file__),'input_sha256':sha(args.input),
        'records':len(records),'effective_records':len(predictions),'row_count':rowcount,'record_scheme_counts':dict(counts),
        'skipped_records':dict(Counter(x['status']for x in records if x['status']!='ready')),
        'frequency_semantics':'Per-position amino-acid reference frequency, organism All, total including gap; not numbering-correctness probability.',
        'frequency_references':frequencies.used,'missing_frequency_positions':len(missing),
        'compatibility_scope':'Exact equality needs independent audit against website references; not a claim of correctness on every possible antibody.'}
    save('manifest.json',manifest)
    print(json.dumps({k:v for k,v in manifest.items()if k!='frequency_references'},ensure_ascii=False))


if __name__=='__main__':main()

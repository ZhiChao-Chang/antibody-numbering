#!/usr/bin/env python3
"""Batch FASTA numbering for antibody variable domains.

Read VH/VL sequences, assign five numbering schemes, and export residue-level
results and run metadata.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import sys

# Resolve defaults from the checkout, so callers may run this script from any
# working directory. Source imports require no package installation.
ROOT=Path(__file__).resolve().parents[1]
if (ROOT/'src').exists():sys.path.insert(0,str(ROOT/'src'))
from antibody_numbering.engine import AnchorNumberer,SCHEMES


def sha(path):
    """Hash an input, software, or reference file for the run manifest."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_fasta(path):
    """Read FASTA records and preserve raw-to-normalized residue correspondence.

    Stable internal record IDs distinguish duplicate headers without deduplicating
    input sequences. Header text is retained verbatim after the leading '>'.
    Malformed FASTA raises an error; invalid sequence records receive skip status.
    """
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
    # source_indices maps each normalized residue back to its one-based offset
    # in the concatenated original sequence, counting '-' characters as positions.
    # Whitespace is removed during parsing and is not a residue coordinate.
    for r in records:
        raw=r['raw_sequence'].upper()
        r['source_indices']=[i+1 for i,aa in enumerate(raw)if aa!='-']
        r['normalized_sequence']=raw.replace('-','').replace('?','X')
        # '*' in a sequence skips that record; a '*' in a FASTA header has no
        # effect. Unknown X and '?' are allowed, but no residue is imputed.
        r['status']=('skipped_stop_character'if'*'in raw else
                     'skipped_invalid_character'if set(raw)-set('ACDEFGHIKLMNPQRSTVWYX?-')else
                     'skipped_empty_sequence'if not r['normalized_sequence']else'ready')
    return records


def main():
    """Run a batch and export residue mappings, summaries, and run metadata."""
    # Stage 1: parse options and load the bundled or supplied framework library.
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--frameworks',type=Path,
        default=ROOT/'data/reference/independent_numbering_public_20261109/framework_templates.json',
        help='Public framework JSON; defaults to the bundled project asset: %(default)s')
    args=p.parse_args()
    # A fresh directory protects existing scientific outputs from accidental reuse.
    assert not args.output.exists(),'Use a new output directory'
    records=read_fasta(args.input);numberer=AnchorNumberer(args.frameworks)
    args.output.mkdir(parents=True)
    columns=['record_id','header','scheme','chain_type','raw_index','normalized_index','input_residue','model_residue','position','status','backend']
    summaries=[];predictions={};counts=Counter();rowcount=0
    # Stage 2: process each record once. The engine returns all five schemes;
    # the CSV expands them into one row per normalized residue per scheme.
    with (args.output/'numbering.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        for record in records:
            if record['status']!='ready':
                # Skipped records remain visible in summaries and the input ledger,
                # even though they produce no residue-level numbering rows.
                for sc in SCHEMES:summaries.append({'record_id':record['record_id'],'scheme':sc,'status':record['status'],'numbered_residues':0})
                continue
            seq=record['normalized_sequence'];result=numberer.number(seq);predictions[record['record_id']]=result
            for sc in SCHEMES:
                value=result['schemes'][sc]
                # The unclassified engine result uses blank lists; classified
                # results carry per-scheme dictionaries and possible failure notes.
                if isinstance(value,list):positions=value;status='unclassified'
                else:positions=value['positions'];status=value['status']
                assert len(positions)==len(seq)
                # Within a record/scheme, every assigned label must be unique.
                # Blank labels are permitted for domain-external or partial output.
                mapped=[x for x in positions if x];assert len(mapped)==len(set(mapped)),(record['record_id'],sc,'duplicate_position')
                counts[sc+':'+status]+=1
                summaries.append({'record_id':record['record_id'],'scheme':sc,'status':status,'numbered_residues':len(mapped)})
                # i is zero-based internally; exported indices below are one-based.
                for i,(aa,position)in enumerate(zip(seq,positions)):
                    writer.writerow({'record_id':record['record_id'],'header':record['header'],'scheme':sc,
                        'chain_type':result.get('chain',''),'raw_index':record['source_indices'][i],'normalized_index':i+1,
                        'input_residue':record['raw_sequence'][record['source_indices'][i]-1],'model_residue':aa,'position':position,
                        'status':status if position or status!='numbered'else'unmapped_position','backend':'own_anchor_region_v1'})
                    rowcount+=1
    # Stage 3: save machine-readable summaries and full intermediate results.
    # The ledger includes the original input sequences and is private run output.
    with (args.output/'records_summary.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=['record_id','scheme','status','numbered_residues']);writer.writeheader();writer.writerows(summaries)
    def save(name,value):(args.output/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    save('input_ledger.json',records);save('backend_outputs.json',predictions)
    import inspect
    engine_path=Path(inspect.getfile(AnchorNumberer))
    # Stage 4: fingerprint the code, input, and framework library.
    manifest={'phase':'COMPLETED','created_at_utc':datetime.now(timezone.utc).isoformat(),
        'backend':'Antibody Numbering',
        'public_frameworks':str(args.frameworks),'public_frameworks_sha256':sha(args.frameworks),
        'engine_source_sha256':sha(engine_path),'runner_source_sha256':sha(__file__),'input_sha256':sha(args.input),
        'records':len(records),'effective_records':len(predictions),'row_count':rowcount,'record_scheme_counts':dict(counts),
        'skipped_records':dict(Counter(x['status']for x in records if x['status']!='ready')),
    }
    save('manifest.json',manifest)
    print(json.dumps(manifest,ensure_ascii=False))


if __name__=='__main__':main()

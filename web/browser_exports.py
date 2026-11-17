"""Browser presentation annotations and per-record ZIP exports.

This module consumes existing CLI outputs. It never changes the numbering
engine, its mappings, or the combined CSV returned by the browser download.
"""
from collections import defaultdict
from pathlib import Path
import base64
import csv
import hashlib
import io
import json
import re
import unicodedata
import zipfile


def csv_text(rows, columns):
    """Serialize filtered rows with the original column order and CSV dialect."""
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def json_text(value):
    return json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n'


def cdr_region(position, chain, definition):
    """Classify labels in the definition's coordinate scheme, including inserts.

    None means that the source scheme did not map this residue. This is distinct
    from a mapped framework residue. Decimal and alphabetic insertion labels
    retain the region of their base position; they are never parsed as floats.
    """
    match = re.fullmatch(r'[HLK]?(\d+)(?:[A-Z]+|\.\d+)?', position or '')
    if not match or chain not in ('H', 'L', 'K'):
        return None
    number = int(match.group(1))
    for index, (start, end) in enumerate(definition['ranges']['H' if chain == 'H' else 'L'], 1):
        if start <= number <= end:
            return 'CDR' + str(index)
    return 'FR'


def record_directory(record):
    """Use stable unique IDs plus a readable, path-safe header fragment."""
    header = unicodedata.normalize('NFC', record['header'])
    fragment = re.sub(r'[^\w.-]+', '_', header)[:60].strip('._-') or 'sequence'
    # CLI IDs are generated internally. Validate instead of trusting header text.
    if not re.fullmatch(r'seq_\d+', record['record_id']):
        raise ValueError('Unexpected record ID')
    return 'sequences/' + record['record_id'] + '__' + fragment


def prepare_browser_result(output_path, definitions_path):
    """Return UI rows and a ZIP with shared metadata and isolated record outputs."""
    output_path = Path(output_path)
    rules_text = Path(definitions_path).read_text(encoding='utf-8')
    rules = json.loads(rules_text)
    files = {p.name: p.read_bytes().decode('utf-8') for p in sorted(output_path.iterdir()) if p.is_file()}
    records = json.loads(files['input_ledger.json'])
    backends = json.loads(files['backend_outputs.json'])
    manifest = json.loads(files['manifest.json'])
    reader = csv.DictReader(io.StringIO(files['numbering.csv']))
    rows = list(reader)
    row_columns = reader.fieldnames
    reader = csv.DictReader(io.StringIO(files['records_summary.csv']))
    summaries = list(reader)
    summary_columns = reader.fieldnames
    grouped_rows, grouped_summaries = defaultdict(list), defaultdict(list)
    source_labels = {}
    for row in rows:
        grouped_rows[row['record_id']].append(row)
        source_labels[(row['record_id'], row['normalized_index'], row['scheme'])] = row
    for row in summaries:
        grouped_summaries[row['record_id']].append(row)

    archive_files = {'manifest.json': files['manifest.json'], 'cdr_definitions.json': rules_text}
    index_rows = []
    for record in records:
        rid = record['record_id']
        directory = record_directory(record)
        local_rows = grouped_rows[rid]
        local_summaries = grouped_summaries[rid]
        archive_files.update({
            directory + '/numbering.csv': csv_text(local_rows, row_columns),
            directory + '/records_summary.csv': csv_text(local_summaries, summary_columns),
            directory + '/backend_outputs.json': json_text({rid: backends[rid]} if rid in backends else {}),
            directory + '/input_ledger.json': json_text([record]),
            directory + '/input.fasta': '>' + record['header'] + '\n' + record['raw_sequence'] + '\n',
        })
        index_rows.append({
            'record_id': rid, 'header': record['header'], 'directory': directory,
            'input_status': record['status'], 'normalized_length': len(record['normalized_sequence']),
            **{row['scheme'] + '_status': row['status'] for row in local_summaries},
        })
    index_columns = ['record_id', 'header', 'directory', 'input_status', 'normalized_length'] + [
        scheme + '_status' for scheme in rules['defaults']]
    archive_files['records_index.csv'] = csv_text(index_rows, index_columns)
    archive_files['README.txt'] = (
        'Antibody Numbering - per-sequence results\n\n'
        'Start with records_index.csv to locate each FASTA record and its processing statuses.\n'
        'Each sequences/ subdirectory contains only that record: numbering.csv (all five\n'
        'schemes), records_summary.csv, backend_outputs.json, input_ledger.json, and input.fasta.\n'
        'Skipped inputs retain their own directory\n'
        'and statuses; their numbering.csv contains the header only. Duplicate FASTA\n'
        'headers remain separate through stable seq_ IDs.\n\n'
        'The outer manifest.json is the original batch provenance, with shared engine,\n'
        'framework and run fingerprints. Its paths refer to the browser virtual filesystem.\n'
        'cdr_definitions.json records the display rules once for the entire batch.\n'
        'export_manifest.json lists export file hashes and the combined CSV fingerprint.\n\n'
        'The separate Download CSV button returns one combined CSV for every sequence\n'
        'and every scheme, unchanged from the CLI. Concatenating the per-record CSV rows\n'
        'in records_index.csv order reproduces those combined data rows.\n'
    )
    archive_files['export_manifest.json'] = json_text({
        'schema': 'per_sequence_zip_v2', 'records': len(records),
        'combined_csv_sha256': hashlib.sha256(files['numbering.csv'].encode()).hexdigest(),
        'exporter_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'files': {name: {'bytes': len(content.encode('utf-8')), 'sha256': hashlib.sha256(content.encode('utf-8')).hexdigest()}
                  for name, content in archive_files.items()},
    })
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as bundle:
        for name, content in archive_files.items():
            bundle.writestr(name, content)

    # Annotate display rows after serialization so exports keep CLI columns intact.
    for row in rows:
        memberships = {}
        for name, definition in rules['definitions'].items():
            source = source_labels.get((row['record_id'], row['normalized_index'], definition['coordinate_scheme']))
            memberships[name] = cdr_region(source['position'], source['chain_type'], definition) if source else None
        row['cdr_regions'] = memberships
    return {'files': files, 'rows': rows, 'summaries': summaries, 'records': records,
            'manifest': manifest, 'cdr_rules': rules,
            'zip': base64.b64encode(archive.getvalue()).decode('ascii')}

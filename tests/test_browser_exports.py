"""Validate browser exports against real CLI output and published CDR policies."""
from pathlib import Path, PurePosixPath
import base64
import csv
import hashlib
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('browser_exports', ROOT / 'web/browser_exports.py')
exports = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exports)
RULES = ROOT / 'web/cdr_definitions.json'


class BrowserExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.temporary.name)
        text = (ROOT / 'examples/6vja_fv.fasta').read_text()
        # Duplicate, non-ASCII and path-like headers must stay distinct and safe.
        text = text.replace('public_6VJA_H', '同名/../../name').replace('public_6VJA_L', '同名/../../name')
        text += '\n>../stop\\name\nACD*EFG\n>empty\n>control\nACDEFGHIKLMNPQRSTVWY\n'
        input_path = cls.directory / 'input.fasta'
        input_path.write_text(text)
        cls.output = cls.directory / 'output'
        subprocess.run([sys.executable, '-B', '-S', str(ROOT / 'scripts/number_antibodies_independent.py'),
                        '--input', str(input_path), '--output', str(cls.output)],
                       check=True, capture_output=True, text=True)
        cls.result = exports.prepare_browser_result(cls.output, RULES)
        cls.archive = zipfile.ZipFile(io.BytesIO(base64.b64decode(cls.result['zip'])))
        cls.index = list(csv.DictReader(io.StringIO(cls.archive.read('records_index.csv').decode())))

    @classmethod
    def tearDownClass(cls):
        cls.archive.close()
        cls.temporary.cleanup()

    def test_combined_csv_is_unchanged(self):
        self.assertEqual(self.result['files']['numbering.csv'].encode(),
                         (self.output / 'numbering.csv').read_bytes())

    def test_split_rows_reassemble_to_the_original_csv(self):
        assembled = []
        for entry in self.index:
            rows = list(csv.DictReader(io.StringIO(self.archive.read(entry['directory'] + '/numbering.csv').decode())))
            self.assertTrue(all(row['record_id'] == entry['record_id'] for row in rows))
            assembled.extend(rows)
        source = csv.DictReader(io.StringIO(self.result['files']['numbering.csv']))
        self.assertEqual(assembled, list(source))
        self.assertEqual(exports.csv_text(assembled, source.fieldnames), self.result['files']['numbering.csv'])

    def test_headers_never_create_unsafe_or_duplicate_paths(self):
        self.assertEqual(len(self.index), 5)
        self.assertEqual(self.index[0]['header'], self.index[1]['header'])
        self.assertEqual(len({entry['directory'] for entry in self.index}), 5)
        for name in self.archive.namelist():
            path = PurePosixPath(name)
            self.assertFalse(path.is_absolute())
            self.assertNotIn('..', path.parts)
            self.assertNotIn('\\', name)

    def test_skipped_inputs_have_their_own_status_and_empty_csv(self):
        for entry in self.index[2:4]:
            directory = entry['directory']
            rows = list(csv.DictReader(io.StringIO(self.archive.read(directory + '/numbering.csv').decode())))
            self.assertEqual(rows, [])
            summary = list(csv.DictReader(io.StringIO(self.archive.read(directory + '/records_summary.csv').decode())))
            self.assertEqual(len(summary), 5)
            self.assertTrue(all(row['status'].startswith('skipped_') for row in summary))
            ledger = json.loads(self.archive.read(directory + '/input_ledger.json'))
            self.assertEqual(ledger[0]['record_id'], entry['record_id'])

    def test_shared_manifest_and_export_hashes(self):
        self.assertEqual(self.archive.read('manifest.json'), (self.output / 'manifest.json').read_bytes())
        self.assertNotIn('numbering.csv', self.archive.namelist())
        manifest = json.loads(self.archive.read('export_manifest.json'))
        self.assertEqual(manifest['schema'], 'per_sequence_zip_v2')
        expected = {'README.txt', 'records_index.csv', 'manifest.json',
                    'export_manifest.json', 'cdr_definitions.json'}
        for record in self.index:
            expected.update(record['directory'] + '/' + name for name in (
                'numbering.csv', 'records_summary.csv', 'backend_outputs.json',
                'input_ledger.json', 'input.fasta',
            ))
        self.assertEqual(set(self.archive.namelist()), expected)
        for name, item in manifest['files'].items():
            content = self.archive.read(name)
            self.assertEqual(len(content), item['bytes'])
            self.assertEqual(hashlib.sha256(content).hexdigest(), item['sha256'])

    def test_cdr_insertions_and_consensus_boundaries(self):
        rules = json.loads(RULES.read_text())['definitions']
        for position, chain, definition, expected in [
            ('H35A', 'H', 'kabat', 'CDR1'), ('H36', 'H', 'kabat', 'FR'),
            ('H112.1', 'H', 'imgt', 'CDR3'), ('L30A', 'L', 'chothia', 'CDR1'),
            ('H95', 'H', 'chothia', 'FR'), ('H96', 'H', 'chothia', 'CDR3'),
            ('L24', 'K', 'chothia', 'FR'), ('L26', 'K', 'chothia', 'CDR1'),
            ('', 'H', 'imgt', None),
        ]:
            self.assertEqual(exports.cdr_region(position, chain, rules[definition]), expected)

    def test_aho_uses_residue_correspondence_not_aho_numeric_ranges(self):
        rows = self.result['rows']
        kabat = next(row for row in rows if row['record_id'] == 'seq_000001' and row['scheme'] == 'kabat' and row['position'] == 'H31')
        aho = next(row for row in rows if row['record_id'] == kabat['record_id'] and row['scheme'] == 'aho' and row['normalized_index'] == kabat['normalized_index'])
        self.assertEqual(aho['cdr_regions']['kabat'], 'CDR1')
        self.assertNotEqual(aho['position'], 'H31')
        # A residue can belong to IMGT CDR1 and Kabat framework at once.
        self.assertTrue(any(row['scheme'] == 'aho' and row['cdr_regions']['imgt'] == 'CDR1'
                            and row['cdr_regions']['kabat'] == 'FR' for row in rows))


if __name__ == '__main__':
    unittest.main()

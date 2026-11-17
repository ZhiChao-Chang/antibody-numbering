from pathlib import Path
import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/number_antibodies_independent.py'
FRAMEWORK = ROOT / 'data/reference/independent_numbering_public_20261109/framework_templates.json'


class CommandLineTests(unittest.TestCase):
    def run_cli(self, output, *extra, cwd=None):
        return subprocess.run(
            [sys.executable, '-B', '-S', str(SCRIPT),
             '--input', str(ROOT / 'examples/6vja_fv.fasta'),
             '--output', str(output), *extra],
            cwd=cwd, capture_output=True, text=True,
        )

    def test_default_framework_resolves_outside_project(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'output'
            result = self.run_cli(output, cwd=directory)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads((output / 'manifest.json').read_text())
            self.assertEqual(Path(manifest['public_frameworks']), FRAMEWORK)
            self.assertEqual(manifest['records'], 2)
            self.assertEqual(manifest['phase'], 'COMPLETED')
            with (output / 'records_summary.csv').open() as stream:
                summaries = list(csv.DictReader(stream))
            self.assertEqual(len(summaries), 10)
            self.assertTrue(all(row['status'] == 'numbered' for row in summaries))
            self.assertTrue(all(int(row['numbered_residues']) > 0 for row in summaries))

    def test_custom_framework_is_used_and_preserves_numbering(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            custom = directory / 'custom_frameworks.json'
            shutil.copyfile(FRAMEWORK, custom)
            output = directory / 'custom_output'
            baseline = directory / 'default_output'
            result = self.run_cli(output, '--frameworks', str(custom), cwd=directory)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads((output / 'manifest.json').read_text())
            self.assertEqual(Path(manifest['public_frameworks']), custom)
            result = self.run_cli(baseline, cwd=directory)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((output / 'numbering.csv').read_bytes(),
                             (baseline / 'numbering.csv').read_bytes())

    def test_numbering_export_columns_and_files(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'output'
            result = self.run_cli(output)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (output / 'numbering.csv').open() as stream:
                reader = csv.DictReader(stream)
                self.assertEqual(reader.fieldnames, [
                    'record_id', 'header', 'scheme', 'chain_type', 'raw_index',
                    'normalized_index', 'input_residue', 'model_residue',
                    'position', 'status', 'backend',
                ])
                self.assertTrue(list(reader))
            self.assertEqual({path.name for path in output.iterdir()}, {
                'numbering.csv', 'records_summary.csv', 'backend_outputs.json',
                'input_ledger.json', 'manifest.json',
            })

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'output'
            output.mkdir()
            sentinel = output / 'keep.txt'
            sentinel.write_text('preserve this file')
            result = self.run_cli(output)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(sentinel.read_text(), 'preserve this file')
            self.assertEqual(list(output.iterdir()), [sentinel])


if __name__ == '__main__':
    unittest.main()

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
            self.assertEqual((output / 'numbering_with_frequencies.csv').read_bytes(),
                             (baseline / 'numbering_with_frequencies.csv').read_bytes())

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

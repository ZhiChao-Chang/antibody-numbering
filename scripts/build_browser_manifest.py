"""Build or verify the browser's versioned, hash-checked Python asset bundle."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RELEASE = '2026.11.17'
PATHS = (
    'src/antibody_numbering/engine.py',
    'src/antibody_numbering/__init__.py',
    'scripts/number_antibodies_independent.py',
    'data/reference/independent_numbering_public_20261109/framework_templates.json',
    'web/browser_exports.py',
    'web/cdr_definitions.json',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify without changing files.')
    args = parser.parse_args()
    expected = {
        'schema': 'browser_algorithm_bundle_v1',
        'release': RELEASE,
        'files': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in PATHS},
    }
    declaration = "const ALGORITHM_RELEASE = '" + RELEASE + "';"
    for name in ('web/app.js', 'web/numbering-worker.mjs'):
        if declaration not in (ROOT / name).read_text():
            raise SystemExit('Browser release mismatch: ' + name)
    if 'web/app.js?v=' + RELEASE not in (ROOT / 'index.html').read_text():
        raise SystemExit('HTML application version does not match the bundle')
    path = ROOT / 'web/algorithm-manifest.json'
    if args.check:
        if not path.exists() or json.loads(path.read_text()) != expected:
            raise SystemExit('Browser asset hashes changed; rebuild algorithm-manifest.json')
        print('Browser bundle verified: ' + RELEASE)
    else:
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(expected, indent=2) + '\n')
        temporary.replace(path)
        print('Browser bundle written: ' + RELEASE)


if __name__ == '__main__':
    main()

# Browser application

[Open the tool](https://zhichao-chang.github.io/antibody-numbering/). Paste VH/VL FASTA, select a sequence and numbering scheme, and download CSV or ZIP. Computation runs in a Web Worker using the same Python CLI as the command-line tool.

## Local preview

From the repository root:

```bash
python3 -m http.server 8000
```

Open `http://localhost:8000/` in a current desktop browser.

## CDR definitions

The numbering scheme and CDR definition can be selected separately. CDR1, CDR2, and CDR3 are blue, orange, and purple.

| Numbering view | Default CDR definition |
| --- | --- |
| IMGT | IMGT |
| Kabat | Kabat |
| Chothia | Chothia consensus (Martin group, 2021) |
| Martin | Chothia consensus, projected by residue index |
| AHo | Kabat, projected by residue index |

IMGT, Kabat, Chothia consensus, and AbM regions are selectable in every view. Exact ranges and source links are in [cdr_definitions.json](cdr_definitions.json).

## Downloads

`numbering.csv` contains every record and all five schemes. ZIP places each record under `sequences/`, with a shared index and run information at the top level. Duplicate headers are kept as separate records, and skipped records retain their input and processing status.

## Deployment and updates

GitHub Pages serves the root of `main`. Keep `.nojekyll`, the Python source, reference files, and `web/runtime/` available at their repository-relative paths.

After modifying bundled Python or JSON assets:

```bash
python3 scripts/build_browser_manifest.py
python3 scripts/build_browser_manifest.py --check
```

For a new release, update the release identifiers in the builder, `app.js`, `numbering-worker.mjs`, and `index.html` together. The worker verifies loaded assets against `algorithm-manifest.json` before running.

## Runtime licenses

The self-hosted runtime is Pyodide 314.0.8. Its source inventory is in [runtime-manifest.json](runtime/runtime-manifest.json). Bundled notices: [Pyodide](runtime/PYODIDE-LICENSE), [Emscripten](runtime/EMSCRIPTEN-LICENSE), and [Python](runtime/PYTHON-LICENSE).

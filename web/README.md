# Browser application

Open https://zhichao-chang.github.io/antibody-numbering/ to paste single- or multi-record FASTA. Each record is one antibody variable domain. All five schemes are computed; the UI initially displays IMGT.

## Execution

`app.js` starts `numbering-worker.mjs` as a module Web Worker. The worker loads the pinned, self-hosted Pyodide 314.0.8 runtime and copies the repository's **unchanged** engine, CLI, package initializer, and default framework JSON into Python's in-memory filesystem. It calls the existing CLI `main()` with `--input` and `--output`; no alternative numbering implementation or remote API is used.

FASTA is passed as a data value, never inserted into Python source. Parsing, numbering, CSV/JSON generation, CDR region projection, and ZIP creation happen inside the worker. `browser_exports.py` consumes the CLI outputs and `cdr_definitions.json` supplies the versioned display rules. The combined CSV remains byte-for-byte unchanged. Network access only retrieves static runtime, algorithm, framework, and optional example files. Input is not stored in browser storage, included in URLs, or submitted to a server. Clearing input terminates the worker; cancellation also terminates it immediately. Repeated runs replace prior in-memory output. Downloaded files remain wherever the user saves them.

The application does not load optional frequency tables. Unmapped, skipped, unsupported, and unclassified statuses follow the CLI. Output files match the CLI; manifest timestamps and virtual filesystem paths naturally differ. The existing engine's scientific validation scope also applies to the browser app.

## Versioned algorithm bundle

The `2026.11.17` application loads the expanded framework library. Versioned URLs and `algorithm-manifest.json` keep the CLI, engine, framework asset, CDR definitions, and export helper synchronized. Before importing the CLI, Python verifies every loaded asset with SHA-256. A partial/stale deployment reports a refresh error instead of silently mixing resources.

After changing any bundled Python/JSON file, run `python3 scripts/build_browser_manifest.py`, then `python3 scripts/build_browser_manifest.py --check`. The runtime version remains pinned separately.

## Local preview

From the repository root:

```bash
python3 -m http.server 8000
```

Visit `http://localhost:8000/`. Use HTTP(S), since module workers and WebAssembly cannot run from a `file://` URL. No build step, JavaScript dependencies, backend, GPU, API key, or login is required. A current desktop browser is recommended for large batches; calculation speed and available memory depend on the client device.

## Deployment

GitHub Pages serves `main` at `/` (repository root). `.nojekyll` preserves the Python source and runtime assets. Keep `src/`, `scripts/`, the default framework path, and `web/runtime/` available at their repository-relative paths. The static files are enough to deploy elsewhere under any base path.

## Runtime attribution

Pyodide 314.0.8 is distributed under the [Mozilla Public License 2.0](runtime/PYODIDE-LICENSE). Runtime sources: https://github.com/pyodide/pyodide/tree/314.0.8 and https://cdn.jsdelivr.net/pyodide/v314.0.8/full/. Exact source URLs and SHA-256 checksums are in `runtime/runtime-manifest.json`. The [Emscripten license](runtime/EMSCRIPTEN-LICENSE) is also included. The [Python license and included component notices](runtime/PYTHON-LICENSE) are provided separately; this repository's MIT license covers its own code, not third-party runtime components.

## CDR definitions

The interface displays the **numbering scheme** and **CDR definition** separately. CDR1, CDR2, and CDR3 use blue, orange, and purple; colors follow the selected definition's residue membership. Region membership is also shown in the data table and tooltips.

Default policies:

| Numbering view | Default CDR definition |
| --- | --- |
| IMGT | IMGT |
| Kabat | Kabat |
| Chothia | Martin group's consensus Chothia definition, updated 28 June 2021 |
| Martin | The same consensus Chothia definition, projected to Martin labels |
| AHo | Kabat definition, projected to AHo labels |

AHo's default is an explicit application policy, not a claim that AHo uniquely specifies Kabat CDRs. Martin numbering and the AbM CDR definition are distinct concepts. Users can independently choose IMGT, Kabat, consensus Chothia, or AbM regions in any numbering view.

The Chothia consensus uses narrower light-chain and H3 boundaries than some older Chothia implementations. The exact inclusive ranges and source URLs are in [cdr_definitions.json](cdr_definitions.json). Sources: [IMGT region definitions](https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html) and the [Martin group's CDR definition tables and consensus notes](http://www.bioinf.org.uk/abs/info.html#cdrdef).

Each definition is evaluated using its own coordinate scheme. Membership is then projected by `(record_id, normalized_index)` to the requested view. Alphabetic and decimal insertions retain their base position's region. Unmapped source positions remain unknown; they are not silently treated as framework. These annotations do not modify the engine's numbering decisions.

## Downloads

**CSV** is one unchanged CLI-compatible file containing every record and all five schemes.

**ZIP** separates sequence-specific files and stores shared run metadata once:

```text
antibody_numbering_results.zip
├── README.txt
├── records_index.csv
├── manifest.json                 # Original batch provenance
├── export_manifest.json          # Export layout and file hashes
├── cdr_definitions.json           # Shared display policies
└── sequences/
    ├── seq_000001__first_header/
    │   ├── numbering.csv          # This record, all five schemes
    │   ├── records_summary.csv
    │   ├── backend_outputs.json
    │   ├── input_ledger.json
    │   ├── input.fasta
    │   └── missing_frequency_positions.json
    └── seq_000002__second_header/
        └── ...
```

Stable record IDs keep duplicate headers separate. Header fragments are sanitized for filesystem paths, while full original names remain in the index and output data. Skipped records keep their own directory, input ledger, five-scheme statuses, and header-only numbering CSV. Concatenating per-record CSV data in index order reproduces the combined CSV.

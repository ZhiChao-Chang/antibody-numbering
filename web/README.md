# Browser application

Open https://zhichao-chang.github.io/antibody-numbering/ to paste single- or multi-record FASTA. Each record is one antibody variable domain. All five schemes are computed; the UI initially displays IMGT.

## Execution

`app.js` starts `numbering-worker.mjs` as a module Web Worker. The worker loads the pinned, self-hosted Pyodide 314.0.8 runtime and copies the repository's **unchanged** engine, CLI, package initializer, and default framework JSON into Python's in-memory filesystem. It calls the existing CLI `main()` with `--input` and `--output`; no alternative numbering implementation or remote API is used.

FASTA is passed as a data value, never inserted into Python source. Parsing, numbering, CSV/JSON generation, and ZIP creation happen inside the worker. Network access only retrieves static runtime, algorithm, framework, and optional example files. Input is not stored in browser storage, included in URLs, or submitted to a server. Clearing input terminates the worker; cancellation also terminates it immediately. Repeated runs replace prior in-memory output. Downloaded files remain wherever the user saves them.

The application does not load optional frequency tables. Unmapped, skipped, unsupported, and unclassified statuses follow the CLI. Output files match the CLI; manifest timestamps and virtual filesystem paths naturally differ. The existing engine's scientific validation scope also applies to the browser app.

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

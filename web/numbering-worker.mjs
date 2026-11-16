// Keep user sequences inside this worker's in-memory Python filesystem.
// Import the repository's unchanged CLI rather than reimplementing its rules.
import { loadPyodide } from './runtime/pyodide.mjs';

let ready;
async function initialize() {
  self.postMessage({type:'progress', message:'正在加载浏览器 Python 环境（首次约 13 MB）…'});
  const py = await loadPyodide({indexURL:new URL('./runtime/', import.meta.url).href});
  const paths = ['src/antibody_numbering/engine.py','src/antibody_numbering/__init__.py',
    'scripts/number_antibodies_independent.py',
    'data/reference/independent_numbering_public_20261109/framework_templates.json'];
  for (const path of paths) {
    const response = await fetch(new URL('../'+path, import.meta.url));
    if (!response.ok) throw new Error('无法加载算法文件：'+path+' ('+response.status+')');
    const bytes = new Uint8Array(await response.arrayBuffer());
    const filename = '/app/'+path;
    py.FS.mkdirTree(filename.slice(0, filename.lastIndexOf('/')));
    py.FS.writeFile(filename, bytes);
  }
  await py.runPythonAsync(`
import sys, importlib.util
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('numbering_cli', '/app/scripts/number_antibodies_independent.py')
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)
`);
  return py;
}

self.onmessage = async ({data}) => {
  const {id, fasta} = data;
  try {
    ready ||= initialize();
    const py = await ready;
    self.postMessage({type:'progress', id, message:'正在编号；所有计算均在此浏览器内进行…'});
    // Passing a value through globals avoids treating input as Python code.
    py.globals.set('fasta_payload', fasta);
    const result = await py.runPythonAsync(`
import io, contextlib, csv, json, zipfile, base64, shutil
from pathlib import Path
input_path = Path('/tmp/input.fasta')
output_path = Path('/tmp/numbering_output')
if output_path.exists(): shutil.rmtree(output_path)
input_path.write_text(fasta_payload, encoding='utf-8')
sys.argv = ['number_antibodies_independent.py', '--input', str(input_path), '--output', str(output_path)]
with contextlib.redirect_stdout(io.StringIO()):
    cli.main()
files = {p.name:p.read_bytes().decode('utf-8') for p in sorted(output_path.iterdir()) if p.is_file()}
archive = io.BytesIO()
with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as bundle:
    for name, content in files.items(): bundle.writestr(name, content)
result = json.dumps({
    'files':files,
    'rows':list(csv.DictReader(io.StringIO(files['numbering_with_frequencies.csv']))),
    'summaries':list(csv.DictReader(io.StringIO(files['records_summary.csv']))),
    'records':json.loads(files['input_ledger.json']),
    'manifest':json.loads(files['manifest.json']),
    'zip':base64.b64encode(archive.getvalue()).decode('ascii')
}, ensure_ascii=False)
# Clear the transient input global; output stays in memory until replaced/closed.
fasta_payload = ''
result
`);
    self.postMessage({type:'result', id, result:JSON.parse(result)});
  } catch (error) {
    self.postMessage({type:'error', id, message:String(error.message || error)});
  }
};

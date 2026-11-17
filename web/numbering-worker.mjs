// Keep user sequences inside this worker's in-memory Python filesystem.
// Import the repository's unchanged CLI rather than reimplementing its rules.
import { loadPyodide } from './runtime/pyodide.mjs';

const ALGORITHM_RELEASE = '2026.11.17.2';
let ready;
async function initialize() {
  self.postMessage({type:'progress', message:'正在加载浏览器 Python 环境（首次约 13 MB）…'});
  const manifestURL = new URL('./algorithm-manifest.json', import.meta.url);
  manifestURL.searchParams.set('release', ALGORITHM_RELEASE);
  const manifestResponse = await fetch(manifestURL, {cache:'no-cache'});
  if (!manifestResponse.ok) throw new Error('算法版本信息加载失败，请刷新页面后重试。');
  const manifest = await manifestResponse.json();
  if (manifest.schema !== 'browser_algorithm_bundle_v1' || manifest.release !== ALGORITHM_RELEASE)
    throw new Error('算法版本正在更新，请刷新页面后重试。');
  const py = await loadPyodide({indexURL:new URL('./runtime/', import.meta.url).href});
  const paths = ['src/antibody_numbering/engine.py','src/antibody_numbering/__init__.py',
    'scripts/number_antibodies_independent.py',
    'data/reference/independent_numbering_public_20261109/framework_templates.json',
    'web/browser_exports.py','web/cdr_definitions.json'];
  for (const path of paths) {
    if (!/^[a-f0-9]{64}$/.test(manifest.files?.[path] || ''))
      throw new Error('算法版本信息不完整，请刷新页面后重试。');
    const resourceURL = new URL('../'+path, import.meta.url);
    resourceURL.searchParams.set('release', ALGORITHM_RELEASE);
    const response = await fetch(resourceURL, {cache:'no-cache'});
    if (!response.ok) throw new Error('无法加载算法文件：'+path+' ('+response.status+')');
    const bytes = new Uint8Array(await response.arrayBuffer());
    const filename = '/app/'+path;
    py.FS.mkdirTree(filename.slice(0, filename.lastIndexOf('/')));
    py.FS.writeFile(filename, bytes);
  }
  py.globals.set('algorithm_bundle_manifest', JSON.stringify(manifest));
  await py.runPythonAsync(`
import sys, importlib.util, hashlib, json
from pathlib import Path
bundle = json.loads(algorithm_bundle_manifest)
for relative, expected in bundle['files'].items():
    actual = hashlib.sha256((Path('/app') / relative).read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError('算法资源版本未同步，请刷新页面后重试：' + relative)
del algorithm_bundle_manifest, bundle
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('numbering_cli', '/app/scripts/number_antibodies_independent.py')
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)
export_spec = importlib.util.spec_from_file_location('browser_exports', '/app/web/browser_exports.py')
exports = importlib.util.module_from_spec(export_spec)
export_spec.loader.exec_module(exports)
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
import io, contextlib, json, shutil
from pathlib import Path
input_path = Path('/tmp/input.fasta')
output_path = Path('/tmp/numbering_output')
if output_path.exists(): shutil.rmtree(output_path)
input_path.write_text(fasta_payload, encoding='utf-8')
sys.argv = ['number_antibodies_independent.py', '--input', str(input_path), '--output', str(output_path)]
with contextlib.redirect_stdout(io.StringIO()):
    cli.main()
result = json.dumps(exports.prepare_browser_result(output_path, '/app/web/cdr_definitions.json'), ensure_ascii=False)
# Clear the transient input global; output stays in memory until replaced/closed.
fasta_payload = ''
result
`);
    self.postMessage({type:'result', id, result:JSON.parse(result)});
  } catch (error) {
    self.postMessage({type:'error', id, message:String(error.message || error)});
  }
};

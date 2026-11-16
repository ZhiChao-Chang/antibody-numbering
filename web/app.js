const $ = id => document.getElementById(id);
let worker, requestId=0, output=null, busy=false;
const statuses={numbered:'已编号',unclassified:'未识别为支持的抗体可变域',unmapped_position:'未映射位置',skipped_stop_character:'跳过：包含终止符 *',skipped_invalid_character:'跳过：包含无效字符',skipped_empty_sequence:'跳过：空序列'};
const describe = value => statuses[value] || value;
const schemeNames={imgt:'IMGT',kabat:'Kabat',chothia:'Chothia',martin:'Martin',aho:'AHo'};
function status(message, kind='') { $('status').textContent=message; $('status').className='status '+kind; }
function showResults(value) { $('results').hidden=!value; if($('guide-number'))$('guide-number').textContent=value?'03':'02'; }
function setBusy(value) { busy=value; $('workspace')?.setAttribute('aria-busy',String(value)); if($('open-file'))$('open-file').disabled=value; $('run').disabled=value; $('cancel').hidden=!value; $('clear').disabled=value; $('example').disabled=value; $('file').disabled=value; $('fasta').readOnly=value; }
function updateCount() { const text=$('fasta').value; const n=(text.match(/^\s*>/gm)||[]).length; $('input-count').textContent=text.trim()?`${n} 条 FASTA 记录 · ${text.length.toLocaleString()} 字符`:'等待输入'; }
function save(name, value, type='text/plain;charset=utf-8') { const link=document.createElement('a'), url=URL.createObjectURL(new Blob([value],{type})); link.href=url; link.download=name; link.click(); setTimeout(()=>URL.revokeObjectURL(url),1500); }
function setupWorker() {
  worker=new Worker(new URL('./numbering-worker.mjs',import.meta.url),{type:'module'});
  worker.onmessage=({data})=>{
    if(data.id!==undefined && data.id!==requestId)return;
    if(data.type==='progress')status(data.message);
    if(data.type==='error') { worker?.terminate(); worker=null; setBusy(false); status('编号未完成：'+data.message,'error'); }
    if(data.type==='result') { output=data.result; setBusy(false); renderResult(); status(`已完成 ${output.records.length} 条记录的处理。结果已在下方显示，可下载全部五种编号体系。`,'success'); }
  };
  worker.onerror=event=>{setBusy(false);status('浏览器计算环境加载失败：'+event.message+'。请检查网络后重试。','error');worker?.terminate();worker=null;};
}
function renderResult() {
  showResults(true); $('record').replaceChildren();
  for(const record of output.records) { const option=document.createElement('option'); option.value=record.record_id; option.textContent=record.record_id+' · '+record.header; $('record').append(option); }
  // Count a record as successful only when every scheme reports success.
  // Review records exclude skipped inputs and include any non-successful scheme.
  const byRecord=new Map();
  for(const row of output.summaries) {
    if(!byRecord.has(row.record_id))byRecord.set(row.record_id,new Map());
    byRecord.get(row.record_id).set(row.scheme,row.status);
  }
  const mapped=output.records.filter(r=>r.status==='ready'&&Object.keys(schemeNames).every(sc=>byRecord.get(r.record_id)?.get(sc)==='numbered')).length;
  const skipped=output.records.filter(x=>x.status!=='ready').length;
  $('summary').replaceChildren();
  for(const [n,text,note] of [[output.records.length,'输入记录','本次输入的全部FASTA记录'],[mapped,'五种体系均已编号','五种体系均返回已编号状态的记录'],[output.records.length-mapped-skipped,'需检查记录','至少一种体系返回未识别、不支持或其他非成功状态；不含跳过记录'],[skipped,'跳过记录','因序列字符或空序列等输入问题而跳过的记录']]) {
    const div=document.createElement('div');div.className='stat';div.title=note;const strong=document.createElement('strong');strong.textContent=n;const span=document.createElement('span');span.textContent=text;div.append(strong,span);$('summary').append(div);
  }
  const counts=$('scheme-statuses');
  if(counts) {
    counts.replaceChildren();
    for(const [scheme,name] of Object.entries(schemeNames)) {
      const item=document.createElement('li'),text=document.createElement('span'),number=document.createElement('strong');
      text.textContent=name;number.textContent=output.summaries.filter(r=>r.scheme===scheme&&r.status==='numbered').length+' / '+output.records.length;
      item.append(text,number);counts.append(item);
    }
  }
  renderRecord();
}
function renderRecord() {
  if(!output)return;
  const id=$('record').value, scheme=$('scheme').value;
  const rows=output.rows.filter(r=>r.record_id===id&&r.scheme===scheme);
  const summary=output.summaries.find(r=>r.record_id===id&&r.scheme===scheme);
  $('record-status').textContent=`${describe(summary?.status||'')} · ${summary?.numbered_residues||0} 个已映射残基`;
  const strip=document.createDocumentFragment(), body=document.createDocumentFragment();
  for(const row of rows) {
    const box=document.createElement('div');box.className='residue';
    const n=Number(row.position.match(/\d+/)?.[0]);
    const isCdr=scheme==='imgt'&&((n>=27&&n<=38)||(n>=56&&n<=65)||(n>=105&&n<=117));
    if(isCdr)box.classList.add('cdr');if(!row.position)box.classList.add('unmapped');
    const aa=document.createElement('b');aa.textContent=row.model_residue;const label=document.createElement('span');label.textContent=row.position||'—';box.append(aa,label);box.title=`原始第${row.raw_index}位；归一化第${row.normalized_index}位；${describe(row.status)}`;strip.append(box);
    const tr=document.createElement('tr');for(const v of [row.raw_index,row.normalized_index,row.model_residue,row.position||'—',describe(row.status)]) {const td=document.createElement('td');td.textContent=v;tr.append(td);}body.append(tr);
  }
  $('residue-strip').replaceChildren(strip);$('rows').replaceChildren(body);
  $('legend').textContent=rows.length?(scheme==='imgt'?'蓝色：IMGT CDR1（27–38）、CDR2（56–65）、CDR3（105–117）；浅色：其余位置。':'上方为当前体系逐残基编号；CDR 着色仅在 IMGT 视图启用。')+' 橙色表示未映射残基。':'此记录没有可展示的逐残基映射，请查看处理状态。';
}
$('fasta').addEventListener('input',updateCount);
$('run').addEventListener('click',()=>{
  const fasta=$('fasta').value;if(!fasta.trim()){status('请先粘贴 FASTA 序列，或载入公开示例。','error');$('fasta').focus();return;}
  if(!fasta.trimStart().startsWith('>')){status('FASTA 每条记录需要以 >名称 开头，请先补上序列名称。','error');return;}
  requestId++;output=null;showResults(false);setBusy(true);status('准备浏览器计算环境…');if(!worker)setupWorker();worker.postMessage({id:requestId,fasta});
});
$('cancel').addEventListener('click',()=>{requestId++;worker?.terminate();worker=null;setBusy(false);status('已取消；输入保留，可重新开始。');});
$('clear').addEventListener('click',()=>{worker?.terminate();worker=null;requestId++;output=null;$('fasta').value='';showResults(false);updateCount();status('已清空输入、结果和浏览器计算内存。');});
$('example').addEventListener('click',async()=>{try{const r=await fetch(new URL('../examples/6vja_fv.fasta',import.meta.url));if(!r.ok)throw Error(r.status);$('fasta').value=await r.text();updateCount();status('已载入 PDB 6VJA 公开的 rituximab VH / VL 两条示例序列。');}catch(e){status('示例加载失败：'+e.message,'error');}});
$('open-file')?.addEventListener('click',()=>$('file').click());
$('file').addEventListener('change',async event=>{const file=event.target.files[0];if(file){$('fasta').value=await file.text();updateCount();status('已在浏览器内读取 '+file.name+'，点击“开始编号”运行。');}event.target.value='';});
$('record').addEventListener('change',renderRecord);$('scheme').addEventListener('change',renderRecord);
$('download-csv').addEventListener('click',()=>{if(output)save('numbering_with_frequencies.csv',output.files['numbering_with_frequencies.csv'],'text/csv;charset=utf-8');});
$('download-zip').addEventListener('click',()=>{if(output)save('antibody_numbering_results.zip',Uint8Array.from(atob(output.zip),c=>c.charCodeAt(0)),'application/zip');});

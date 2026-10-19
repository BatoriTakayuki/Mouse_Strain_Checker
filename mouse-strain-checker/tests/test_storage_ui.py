"""Exercise the actual UI event flow with a deterministic DOM and API fixture."""
from pathlib import Path
import shutil
import subprocess
import unittest


STATIC = Path(__file__).resolve().parents[1] / "mouse_strain" / "static"


JS_CHECKS = r"""
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const ui = require(process.argv[1]);
const plan = {token:'a'.repeat(64),files:3};
assert.equal(ui.canConfirmCleanup(plan,true,false),true);
for (const [p,agreed,busy] of [[plan,false,false],[plan,'true',false],[plan,true,true],
  [null,true,false],[{...plan,files:0},true,false],[{...plan,token:'bad'},true,false]]) {
  assert.equal(ui.canConfirmCleanup(p,agreed,busy),false);
}
assert.ok(ui.sampleExplanation({cleanup:{scope:'cache'}}).includes('中間ファイルは削除済み'));

(async () => {
  class Node {
    constructor(tag='div') { this.tagName=tag; this.childNodes=[]; this.listeners={};
      this.disabled=false; this.checked=false; this.open=false; this.hidden=false;
      this.classList={toggle(){}}; this.style={setProperty(){}}; this.textContent=''; }
    append(...nodes) { this.childNodes.push(...nodes); }
    replaceChildren(...nodes) { this.childNodes=nodes; }
    setAttribute() {} removeAttribute() {}
    addEventListener(name,fn) { (this.listeners[name] ||= []).push(fn); }
    async fire(name) { const event={preventDefault(){}};
      for (const fn of this.listeners[name] || []) await fn(event); }
    showModal() { this.open=true; }
    close() { this.open=false; this.fire('close'); }
  }
  const nodes = new Map();
  const get = id => { if (!nodes.has(id)) nodes.set(id,new Node()); return nodes.get(id); };
  const doc = {getElementById:get,createElement:tag=>new Node(tag),createDocumentFragment:()=>new Node('fragment')};
  const complete = {id:'a'.repeat(32),status:'completed',request:{srx:['SRX1']},created_at:'2026-10-19'};
  const running = {id:'b'.repeat(32),status:'running',request:{srx:['SRX2']},created_at:'2026-10-19'};
  let jobs=[complete,running], cacheBytes=4096, cleanupCalls=0, pendingCleanup=null;
  const calls=[];
  const fetch = async (path,options={}) => {
    const body=options.body ? JSON.parse(options.body) : null; calls.push({path,body});
    let data;
    if (path==='/api/health') data={tools:{},reference_ready:true};
    else if (path==='/api/jobs') data=jobs;
    else if (path==='/api/storage') data={free_bytes:1024,total_bytes:8192,
      jobs:jobs.map(j=>({id:j.id,status:j.status,accessions:j.request.srx,total_bytes:8192,
        cache_bytes:j===complete ? cacheBytes : null,deletable:j===complete})),reference:{deletable:false}};
    else if (path==='/api/storage/preview') data={...plan,scope:body.scope,job_id:body.job_id,
      title:'SRX1',effect:'結果を保持',bytes:cacheBytes,paths:['synthetic/cache/SRR1.bam']};
    else if (path==='/api/storage/cleanup') {
      cleanupCalls++; await new Promise(resolve=>{pendingCleanup=resolve;});
      cacheBytes=0; data={bytes:4096,message:'対象ファイルを削除しました。'};
    } else throw new Error('Unexpected API path: '+path);
    return {ok:true,text:async()=>JSON.stringify(data)};
  };
  const window={addEventListener(){},setInterval(){},setTimeout,clearTimeout};
  vm.runInNewContext(fs.readFileSync(process.argv[1],'utf8'),
    {document:doc,window,location:{hash:'',pathname:'/',search:''},history:{replaceState(){}},
      fetch,AbortController,console,Date,Error,TypeError});
  const flush=async()=>{for(let i=0;i<8;i++) await new Promise(setImmediate);};
  await flush(); get('storage-panel').open=true; await get('storage-panel').fire('toggle'); await flush();
  const rows=get('storage-table').childNodes[1].childNodes;
  const inactiveActions=rows[0].childNodes[4].childNodes;
  const activeActions=rows[1].childNodes[4].childNodes;
  assert.equal(inactiveActions[0].disabled,false);
  assert.equal(activeActions[0].disabled,true); assert.equal(activeActions[1].disabled,true);
  await inactiveActions[0].fire('click'); await flush();
  assert.equal(get('cleanup-dialog').open,true); assert.equal(get('confirm-cleanup').disabled,true);
  await get('cleanup-form').fire('submit'); assert.equal(cleanupCalls,0);
  await get('cancel-cleanup').fire('click'); assert.equal(get('cleanup-dialog').open,false);
  assert.equal(cleanupCalls,0);
  await inactiveActions[0].fire('click'); await flush();
  get('cleanup-agree').checked=true; await get('cleanup-agree').fire('change');
  assert.equal(get('confirm-cleanup').disabled,false);
  const deletion=get('cleanup-form').fire('submit'); await flush();
  assert.equal(get('confirm-cleanup').disabled,true); assert.equal(get('cancel-cleanup').disabled,true);
  await get('cleanup-form').fire('submit'); assert.equal(cleanupCalls,1);
  assert.deepEqual(calls.find(call=>call.path==='/api/storage/cleanup').body,{token:plan.token,confirm:true});
  pendingCleanup(); await deletion; await flush();
  assert.equal(get('cleanup-dialog').open,false);
  assert.ok(get('storage-feedback').textContent.includes('4 KiB'));
  assert.equal(get('storage-table').childNodes[1].childNodes[0].childNodes[4].childNodes[0].disabled,true);
  console.log('Storage confirmation, cancel, active-job protection and duplicate-submit checks passed.');
})().catch(error=>{console.error(error);process.exitCode=1;});
"""


class StorageUITests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("node"), "Node is optional; JS_CHECKS can run in an available JS runtime")
    def test_confirmation_and_delete_event_flow(self):
        result = subprocess.run([shutil.which("node"), "-e", JS_CHECKS, str(STATIC / "app.js")],
                                text=True, encoding="utf-8", capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()

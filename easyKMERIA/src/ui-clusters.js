/* ===== 集群连接 ===== */
const SCHED_NAME = {slurm:'SLURM', sge:'SGE', pbs:'PBS', local:tr('直接运行')};
let editingId = null, testResults = {};

async function loadClusters(){
  $('clustersStatic').hidden = state.backend;
  $('btnAddCluster').disabled = !state.backend;
  if (state.backend) {
    try { state.clusters = (await api('GET', '/api/clusters')).clusters; } catch (e) { toast(e.message); }
  }
  renderClusters();
}
function renderClusters(){
  const box = $('clusterList');
  if (!state.backend) { box.innerHTML = ''; return; }
  if (!state.clusters.length) { box.innerHTML = tr('<div class="empty">还没有保存的连接。点右上角“添加连接”，填上服务器地址和账号。</div>'); return; }
  box.innerHTML = state.clusters.map(c => {
    const authTxt = c.auth === 'password' ? (tr('密码 · ') + (c.has_password ? tr('已保存') : tr('未保存'))) : (tr('私钥 · ') + (c.has_key ? tr('已保存') : tr('未保存')));
    const r = testResults[c.id];
    return '<article class="card cl-card" data-id="' + esc(c.id) + '">' +
      '<div class="cl-top"><h2>' + esc(c.name) + '</h2><span>' + (c.id === state.activeId ? tr('<span class="chip ok">当前</span> ') : '') + '<span class="chip">' + esc(SCHED_NAME[c.scheduler] || c.scheduler) + '</span></span></div>' +
      '<p class="mono muted">' + esc(c.username) + '@' + esc(c.host) + ':' + c.port + '</p>' +
      tr('<div class="kv"><div><span class="lab">工作路径</span><span>') + esc(c.workdir || tr('未设置')) + tr('</span></div><div><span class="lab">队列</span><span>') + esc(c.queue || tr('未设置')) + '</span></div>' +
      tr('<div><span class="lab">认证</span><span>') + esc(authTxt) + tr('</span></div><div><span class="lab">安装目录</span><span>') + esc(c.kmeria_home || tr('未设置')) + '</span></div></div>' +
      (c.host_fp ? tr('<p class="hint mono">服务器指纹 ') + esc(c.host_fp) + '</p>' : tr('<p class="hint">还没有记录服务器指纹，第一次测试连接时会自动记录。</p>')) +
      '<div class="cl-res" ' + (r ? '' : 'hidden') + ' id="res_' + esc(c.id) + '">' + (r || '') + '</div>' +
      tr('<div class="cl-actions"><button type="button" class="btn" data-a="test">测试连接</button><button type="button" class="btn" data-a="use">设为当前</button>') +
      tr('<button type="button" class="btn" data-a="edit">编辑</button><button type="button" class="btn danger" data-a="del">删除</button></div></article>');
  }).join('');
}
function testHtml(r){
  const row = (lvl, t) => '<div class="chk"><span class="pill ' + lvl + '">' + ({ok:tr('正常'), warn:tr('注意'), err:tr('错误')})[lvl] + '</span><span>' + esc(t) + '</span></div>';
  let h = row('ok', tr('登录成功：') + r.login);
  h += r.fingerprint_new ? row('ok', tr('已记录服务器指纹 ') + r.fingerprint + tr('，以后指纹变化会拒绝连接')) : row('ok', tr('服务器指纹与上次一致'));
  const sf = {slurm:tr('SLURM（sbatch）'), qsub:tr('SGE 或 PBS（qsub）'), none:tr('没有检测到调度系统')}[r.scheduler_found];
  h += row(r.scheduler_found === 'none' ? 'warn' : 'ok', tr('调度系统：') + sf);
  if (r.queues && r.queues.length) h += row(r.queue_exists === false ? 'warn' : 'ok', tr('可用队列：') + r.queues.join('、') + (r.queue_exists === false ? tr('（你填的队列不在其中）') : ''));
  if ('workdir_exists' in r) h += row(r.workdir_exists ? 'ok' : 'warn', r.workdir_exists ? tr('工作路径存在') : tr('工作路径不存在，部署时会自动创建'));
  if ('kmeria_found' in r) h += row(r.kmeria_found ? 'ok' : 'warn', r.kmeria_found ? tr('找到了 KMERIA（kmeria_wrapper.pl）') : tr('安装目录里没有找到 kmeria_wrapper.pl，请到“流程生成”的“安装 KMERIA”页按步骤安装'));
  return h;
}

function openClusterForm(c){
  editingId = c ? c.id : null;
  $('cmTitle').textContent = c ? tr('编辑连接') : tr('添加连接');
  $('cf_name').value = c ? c.name : ''; $('cf_host').value = c ? c.host : ''; $('cf_port').value = c ? c.port : 22;
  $('cf_user').value = c ? c.username : ''; $('cf_auth').value = c ? c.auth : 'password'; $('cf_scheduler').value = c ? c.scheduler : 'slurm';
  $('cf_workdir').value = c ? c.workdir : ''; $('cf_queue').value = c ? c.queue : ''; $('cf_home').value = c ? c.kmeria_home : '';
  $('cf_password').value = ''; $('cf_key').value = ''; $('cf_pass').value = ''; $('cf_reset').checked = false;
  $('cf_pwHint').textContent = c && c.has_password ? tr('已保存密码，留空表示不修改。') : '';
  $('cf_keyHint').textContent = c && c.has_key ? tr('已保存私钥，留空表示不修改。') : '';
  $('cf_resetBox').hidden = !c; $('cf_err').hidden = true;
  syncAuth(); openModal('clusterModal');
}
function syncAuth(){ const k = $('cf_auth').value === 'key'; $('cf_pwBox').hidden = k; $('cf_keyBox').hidden = !k; }

async function submitCluster(){
  const body = {name:$('cf_name').value.trim(), host:$('cf_host').value.trim(), port:Number($('cf_port').value) || 22, username:$('cf_user').value.trim(),
    auth:$('cf_auth').value, scheduler:$('cf_scheduler').value, queue:$('cf_queue').value.trim(), workdir:$('cf_workdir').value.trim(), kmeria_home:$('cf_home').value.trim(),
    password:$('cf_password').value, key_text:$('cf_key').value, passphrase:$('cf_pass').value, reset_fingerprint:$('cf_reset').checked};
  const btn = $('cf_save'); btn.disabled = true; $('cf_err').hidden = true;
  try {
    const r = editingId ? await api('PUT', '/api/clusters/' + editingId, body) : await api('POST', '/api/clusters', body);
    if (!state.activeId) { state.activeId = r.cluster.id; LS.set('easykmera.active', state.activeId); }
    closeModal('clusterModal'); toast(tr('已保存')); await loadClusters();
  } catch (e) { $('cf_err').textContent = e.message; $('cf_err').hidden = false; }
  finally { btn.disabled = false; }
}

function initClusters(){
  $('btnAddCluster').addEventListener('click', () => openClusterForm(null));
  $('cf_auth').addEventListener('change', syncAuth);
  $('clusterForm').addEventListener('submit', submitCluster);
  $('clusterList').addEventListener('click', async e => {
    const b = e.target.closest('[data-a]'); if (!b) return;
    const card = b.closest('[data-id]'), id = card.dataset.id, c = state.clusters.find(x => x.id === id);
    if (b.dataset.a === 'edit') return openClusterForm(c);
    if (b.dataset.a === 'use') { state.activeId = id; LS.set('easykmera.active', id); renderClusters(); toast(tr('已设为当前连接')); return; }
    if (b.dataset.a === 'del') {
      if (b.dataset.armed !== '1') { b.dataset.armed = '1'; b.textContent = tr('确认删除'); b.classList.add('armed'); setTimeout(() => { b.dataset.armed = ''; b.textContent = tr('删除'); b.classList.remove('armed'); }, 4000); return; }
      try { await api('DELETE', '/api/clusters/' + id); if (state.activeId === id) state.activeId = ''; delete testResults[id]; toast(tr('已删除')); await loadClusters(); } catch (er) { toast(er.message); }
      return;
    }
    if (b.dataset.a === 'test') {
      const box = $('res_' + id); box.hidden = false; box.textContent = tr('正在连接…'); b.disabled = true;
      try { const r = await api('POST', '/api/clusters/' + id + '/test', {}); testResults[id] = testHtml(r); box.innerHTML = testResults[id]; if (r.fingerprint_new) await loadClustersQuiet(); }
      catch (er) { testResults[id] = tr('<div class="chk"><span class="pill err">错误</span><span>') + esc(er.message) + '</span></div>'; box.innerHTML = testResults[id]; }
      finally { b.disabled = false; }
    }
  });
  routeHooks.clusters = loadClusters;
}
async function loadClustersQuiet(){ try { state.clusters = (await api('GET', '/api/clusters')).clusters; renderClusters(); } catch (e) {} }

/* ===== AI 设置 ===== */
async function loadAi(){
  $('aiStatic').hidden = state.backend; $('btnAiEdit').disabled = $('btnAiTest').disabled = !state.backend;
  $('dataDir').textContent = state.backend ? state.dataDir : tr('未连接到程序（只读预览）');
  if (state.backend) { try { Object.assign(state.ai, await api('GET', '/api/ai/settings')); } catch (e) { toast(e.message); } }
  renderAi();
}
function renderAi(){
  const a = state.ai;
  $('aiChip').textContent = a.has_key ? tr('已配置') : tr('未设置'); $('aiChip').className = 'chip' + (a.has_key ? ' ok' : '');
  $('aiUrl').textContent = a.base_url || 'https://api.deepseek.com'; $('aiModel').textContent = a.model || 'deepseek-chat';
  $('aiKeyState').textContent = a.has_key ? tr('已保存') : tr('未保存');
}
function openAiForm(){
  $('ai_key').value = ''; $('ai_base').value = state.ai.base_url || 'https://api.deepseek.com'; $('ai_model').value = state.ai.model || 'deepseek-chat';
  $('ai_keyHint').textContent = state.ai.has_key ? tr('已保存 API Key；留空表示不修改。') : tr('到 platform.deepseek.com 创建 API Key 后粘贴到这里。');
  $('ai_key').placeholder = state.ai.has_key ? tr('已保存，留空不修改') : 'sk-…';
  $('ai_ca').value = state.ai.ca_bundle || ''; $('ai_clear').hidden = !state.ai.has_key; $('ai_err').hidden = true; openModal('aiModal');
}
function initAi(){
  $('btnAiEdit').addEventListener('click', openAiForm);
  $('aiForm').addEventListener('submit', async () => {
    $('ai_err').hidden = true; $('ai_save').disabled = true;
    try { Object.assign(state.ai, await api('POST', '/api/ai/settings', {api_key:$('ai_key').value, base_url:$('ai_base').value, model:$('ai_model').value, ca_bundle:$('ai_ca').value})); closeModal('aiModal'); renderAi(); toast(tr('已保存')); }
    catch (e) { $('ai_err').textContent = e.message; $('ai_err').hidden = false; } finally { $('ai_save').disabled = false; }
  });
  $('ai_clear').addEventListener('click', async () => {
    try { Object.assign(state.ai, await api('POST', '/api/ai/settings', {clear_key:true, base_url:$('ai_base').value, model:$('ai_model').value, ca_bundle:$('ai_ca').value})); closeModal('aiModal'); renderAi(); toast(tr('已清除 API Key')); }
    catch (e) { $('ai_err').textContent = e.message; $('ai_err').hidden = false; }
  });
  $('btnAiTest').addEventListener('click', async () => {
    $('aiTestMsg').textContent = tr('正在测试…'); $('btnAiTest').disabled = true;
    try { await api('POST', '/api/ai/test', {}); $('aiTestMsg').textContent = tr('连接正常，DeepSeek 已响应。'); }
    catch (e) { $('aiTestMsg').textContent = e.message; } finally { $('btnAiTest').disabled = false; }
  });
  routeHooks.ai = loadAi;
}

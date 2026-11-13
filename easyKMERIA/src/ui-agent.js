/* ===== AI 助手: 只调用已有的白名单动作; AI 只产出经后端校验的表单字段和文字建议 ===== */
const ag = {busy:false, stop:false, wake:null, pollMs:30000};
const AG_LABEL = {project:tr('项目名'), workDir:tr('工作目录'), kmeriaHome:tr('KMERIA 安装目录'), fastqDir:tr('原始数据目录'), condaEnv:tr('conda 环境'), dataMode:tr('数据来源方式'),
  defaultDepth:tr('默认深度'), readLen:tr('读长'), ploidy:tr('倍性'), kmer:tr('k 值'), assoTool:tr('关联工具'), counter:tr('计数工具'), step:tr('运行步骤'), scheduler:tr('调度系统'),
  queue:tr('队列'), threads:tr('线程'), memory:tr('内存(GB)'), time:tr('时限'), sampleText:tr('样本表'), trait:tr('性状序号')};
const AG_PATH_OK = /^\/[A-Za-z0-9_.\/@+=-]+$/;
const agCluster = () => state.clusters.find(c => c.id === $('ag_sel').value);
const agSleep = ms => new Promise(res => { const t = setTimeout(res, ms); ag.wake = () => { clearTimeout(t); res(); }; });

function agItem(cls, html){
  const e = $('ag_empty'); if (e) e.remove();
  const d = document.createElement('div'); d.className = 'ag-item ' + cls; d.innerHTML = html;
  $('ag_log').appendChild(d); d.scrollIntoView({block:'nearest'}); return d;
}
const agNote = (title, html) => agItem('', '<div class="ag-h"><b>' + esc(title) + '</b></div>' + (html ? '<p>' + html + '</p>' : ''));
const agErr = (title, msg) => agItem('err', tr('<div class="ag-h"><span class="pill err">出错</span><b>') + esc(title) + '</b></div><p>' + esc(msg) + '</p>');
function agStep(title){
  const d = agItem('', tr('<div class="ag-h"><span class="pill warn">进行中</span><b>') + esc(title) + '</b></div><pre class="rlog" hidden></pre>');
  const pill = d.querySelector('.pill'), pre = d.querySelector('pre');
  return {
    log(t){ pre.hidden = false; pre.textContent = t; pre.scrollTop = pre.scrollHeight; },
    ok(t){ pill.className = 'pill ok'; pill.textContent = tr('完成'); if (t != null) this.log(t); },
    fail(t){ pill.className = 'pill err'; pill.textContent = tr('失败'); d.classList.add('err'); if (t != null) this.log(t); },
    skip(t){ pill.className = 'pill'; pill.textContent = tr('已跳过'); if (t != null) this.log(t); }
  };
}
/* 批准卡片: 返回 true(批准) / false(拒绝或停止) */
function agAsk(title, html, okLabel){
  return new Promise(res => {
    const d = agItem('ask', tr('<div class="ag-h"><span class="pill warn">等你批准</span><b>') + esc(title) + '</b></div>' + html +
      '<div class="ag-btns"><button type="button" class="btn primary" data-a="1">' + esc(okLabel || tr('批准')) + tr('</button><button type="button" class="btn" data-a="0">不做了</button></div>'));
    const done = v => { d.querySelectorAll('button').forEach(b => { b.disabled = true; }); d.classList.remove('ask'); d.querySelector('.pill').className = 'pill ' + (v ? 'ok' : ''); d.querySelector('.pill').textContent = v ? tr('已批准') : tr('已拒绝'); ag.wake = null; res(v); };
    d.addEventListener('click', e => { const b = e.target.closest('[data-a]'); if (b) done(b.dataset.a === '1'); });
    ag.wake = () => done(false);
  });
}
function agBegin(){
  if (!state.backend) { toast(tr('需要运行 easyKMERIA 程序')); return null; }
  if (ag.busy) { toast(tr('已有任务在进行，请先等它结束或点“停止”')); return null; }
  ag.busy = true; ag.stop = false;
  $('ag_stop').hidden = false; ['ag_fill','ag_install','ag_diag','ag_auto'].forEach(i => { $(i).disabled = true; });
  return true;
}
function agEnd(){
  ag.busy = false; ag.wake = null; $('ag_stop').hidden = true;
  ['ag_fill','ag_install','ag_diag','ag_auto'].forEach(i => { $(i).disabled = false; });
}
const agNeedAi = () => { if (!state.ai.has_key) { agErr(tr('还没有配置 DeepSeek'), tr('请先到“AI 设置”保存 API Key。')); return false; } return true; };
function agNeedCluster(){
  const c = agCluster(); if (!c) { agErr(tr('没有服务器连接'), tr('请先到“集群连接”添加并测试一个连接。')); return null; }
  if (!c.workdir) { agErr(tr('连接里没有工作路径'), tr('请编辑这个连接，填写工作路径（脚本只会在它之内读写）。')); return null; }
  return c;
}
const agCfgForAi = () => { const c = readCfg(); delete c.sampleText; delete c.phenoText; delete c.preCmd; return c; };
function agCheckNote(){
  const errs = lastGen.checks.filter(c => c.level === 'err');
  if (errs.length) agItem('err', tr('<div class="ag-h"><span class="pill err">还差 ') + errs.length + tr(' 项</span><b>表单检查未通过</b></div><ul>') + errs.slice(0, 6).map(c => '<li>' + esc(c.msg) + '</li>').join('') + tr('</ul><p><a href="#generate">到“流程生成”页补充</a></p>'));
  else agNote(tr('表单检查通过'), tr('可以继续部署了。'));
}

/* ---- ① 一句话填表 ---- */
async function agFill(){
  const text = $('ag_text').value.trim();
  if (!text) { toast(tr('先写一句话描述')); return; }
  if (!agBegin()) return;
  try {
    if (!agNeedAi()) return;
    const st = agStep(tr('让 DeepSeek 把描述转成表单字段'));
    let r; try { r = await api('POST', '/api/agent/fill', {text, current: agCfgForAi()}); } catch (e) { st.fail(e.message); return; }
    st.ok();
    const keys = Object.keys(r.patch), cur = readCfg();
    if (!keys.length) { agNote(tr('没有可以填写的内容'), esc(r.notes || '') + (r.missing.length ? tr('还缺：') + esc(r.missing.join('、')) : '')); return; }
    const rows = keys.map(k => '<tr><td>' + esc(AG_LABEL[k] || k) + '</td><td class="mono">' + esc(k === 'sampleText' ? tr('（') + r.patch[k].split('\n').filter(Boolean).length + tr(' 行）') : String(cur[k] || '—')) + '</td><td class="mono">' + esc(k === 'sampleText' ? tr('（') + r.patch[k].split('\n').filter(Boolean).length + tr(' 行）') : r.patch[k]) + '</td></tr>').join('');
    const ok = await agAsk(tr('应用到“流程生成”表单？'), (r.notes ? '<p>' + esc(r.notes) + '</p>' : '') + tr('<table class="ag-diff"><tr><th>字段</th><th>现在</th><th>改成</th></tr>') + rows + '</table>' +
      (r.missing.length ? tr('<p>还缺：') + esc(r.missing.join('、')) + tr('。请到表单里补充，或再描述一次。</p>') : ''), tr('应用'));
    if (!ok) return;
    loadCfg(Object.assign(readCfg(), r.patch)); $('demo').hidden = true;
    agNote(tr('已填入表单'), tr('改了 ') + keys.length + tr(' 项。')); agCheckNote();
  } finally { agEnd(); }
}

/* ---- ② 检查并安装 KMERIA ---- */
function agInstallScript(home){
  return ['#!/usr/bin/env bash', tr('# easyKMERIA 生成: 安装 KMERIA (官方 README 的预编译方式)。不会使用 sudo, 也不会删除任何文件。'),
    'set -uo pipefail', 'HOME_DIR="' + home + '"',
    'finish() { echo "INSTALL_DONE $1"; exit "$1"; }',
    tr('command -v git >/dev/null 2>&1 || { echo "[错误] 这台机器没有 git"; finish 2; }'),
    'CONDA=""; for t in mamba micromamba conda; do if command -v "$t" >/dev/null 2>&1; then CONDA="$t"; break; fi; done',
    tr('[ -n "$CONDA" ] || { echo "[错误] 没找到 conda / mamba / micromamba, 请先安装 Miniconda 或联系管理员"; finish 2; }'),
    tr('if [ -e "$HOME_DIR" ] && [ -n "$(ls -A "$HOME_DIR" 2>/dev/null)" ]; then echo "[错误] $HOME_DIR 已存在且不是空目录, 为避免覆盖已停止"; finish 3; fi'),
    tr('mkdir -p "$(dirname "$HOME_DIR")" || { echo "[错误] 无法创建上级目录, 请检查权限"; finish 4; }'),
    tr('echo "[1/4] 下载 KMERIA"'), tr('git clone https://github.com/Sh1ne111/KMERIA.git "$HOME_DIR" || { echo "[错误] git clone 失败, 请检查这台机器能否访问 github.com"; finish 5; }'),
    tr('echo "[2/4] 创建 conda 环境 kmeriaenv (需要几分钟)"'), 'cd "$HOME_DIR" || finish 5',
    tr('if "$CONDA" env list 2>/dev/null | grep -Fwq kmeriaenv; then echo "环境 kmeriaenv 已存在, 跳过"; else "$CONDA" env create -f kmeria_env.yaml || { echo "[错误] 创建 conda 环境失败"; finish 6; }; fi'),
    tr('echo "[3/4] 设置库路径和执行权限"'),
    'chmod 755 "$HOME_DIR"/bin/* "$HOME_DIR"/external_tools/* "$HOME_DIR"/bimbamAsso/* 2>/dev/null || true',
    'export LD_LIBRARY_PATH="$HOME_DIR/lib:${LD_LIBRARY_PATH:-}"', 'export PATH="$HOME_DIR/bin:$HOME_DIR/bimbamAsso:$HOME_DIR/external_tools:$PATH"',
    tr('echo "[4/4] 验证"'),
    tr('if command -v kmeria >/dev/null 2>&1; then kmeria 2>&1 | head -n 12; echo INSTALL_OK; finish 0; else echo "[警告] 装完后仍找不到 kmeria 命令, 请核对安装目录里 bin/ 下的内容"; finish 7; fi')].join('\n') + '\n';
}
function agCheckScript(home){
  return ['#!/usr/bin/env bash', 'HOME_DIR="' + home + '"',
    'if [ -d "$HOME_DIR/bin" ]; then echo KMERIA_PRESENT; else echo KMERIA_MISSING; fi',
    'for t in git conda mamba micromamba make gcc; do if command -v "$t" >/dev/null 2>&1; then echo "HAVE $t"; else echo "MISSING $t"; fi; done'].join('\n') + '\n';
}
/* 返回 true 表示 KMERIA 已就绪 */
async function agEnsureKmeria(c){
  const home = readCfg().kmeriaHome.replace(/\/+$/, '');
  if (!AG_PATH_OK.test(home)) { agErr(tr('KMERIA 安装目录不合法'), tr('请在表单第 1 步填写绝对路径（只含字母、数字和 _ . / @ + = -）。')); return false; }
  const dir = joinPath(c.workdir, 'easykmera_setup');
  let st = agStep(tr('检查服务器上有没有 KMERIA'));
  let out;
  try {
    await api('POST', '/api/clusters/' + c.id + '/upload', {dir, files:[{name:'check_env.sh', content:agCheckScript(home)}]});
    out = (await api('POST', '/api/clusters/' + c.id + '/exec', {dir, action:'bash', script:'check_env.sh'})).output;
  } catch (e) { st.fail(e.message); return false; }
  if (/KMERIA_PRESENT/.test(out)) { st.ok(out.trim() + tr('\n\n已经安装，无需再装。')); return true; }
  st.ok(out.trim());
  if (!(/HAVE (conda|mamba|micromamba)/.test(out) && /HAVE git/.test(out))) { agErr(tr('缺少安装条件'), tr('这台服务器上没找到 git 或 conda/mamba，我没法替你装。请联系管理员或手动安装后再来。')); return false; }
  const script = agInstallScript(home);
  if (!await agAsk(tr('安装 KMERIA 到 ') + home, tr('<p>会在登录节点上后台运行下面的脚本：克隆 GitHub 仓库、创建 conda 环境 kmeriaenv。目录已存在且非空时会自动停止，不会覆盖。需要这台服务器能访问 github.com 和 conda 源，耗时几分钟到几十分钟。</p><pre class="rlog">') + esc(script) + '</pre>', tr('批准安装'))) return false;
  st = agStep(tr('安装 KMERIA（后台运行）'));
  try {
    await api('POST', '/api/clusters/' + c.id + '/upload', {dir, files:[{name:'install_kmeria.sh', content:script}]});
    await api('POST', '/api/clusters/' + c.id + '/exec', {dir, action:'bash_bg', script:'install_kmeria.sh'});
  } catch (e) { st.fail(e.message); return false; }
  const t0 = Date.now(); let log = '';
  while (Date.now() - t0 < 45 * 60 * 1000) {
    await agSleep(Math.min(5000, ag.pollMs)); if (ag.stop) { st.skip(tr('已停止监控（服务器上的安装脚本仍在运行，日志在 ') + dir + tr('/install_kmeria.sh.log）')); return false; }
    try { log = (await api('POST', '/api/clusters/' + c.id + '/tail', {path: dir + '/install_kmeria.sh.log', lines: 60})).output; } catch (e) { continue; }
    st.log(log);
    const m = log.match(/INSTALL_DONE (\d+)/); if (m) { if (m[1] === '0') { st.ok(); agNote(tr('KMERIA 安装完成'), tr('已在 ') + esc(home) + tr(' 装好，部署脚本会自动加载它。')); return true; } break; }
  }
  st.fail();
  await agDiagnose(c, 'install', log || tr('安装脚本超过 45 分钟没有结束'));
  return false;
}
async function agInstall(){
  if (!agBegin()) return;
  try { const c = agNeedCluster(); if (c) await agEnsureKmeria(c); } finally { agEnd(); }
}

/* ---- ③ 诊断 ---- */
function agMergePatch(fixes){ const p = {}; fixes.forEach(f => Object.keys(f.patch || {}).forEach(k => { if (!(k in p)) p[k] = f.patch[k]; })); return p; }
function agDiffRows(patch){ const cur = readCfg(); return Object.keys(patch).map(k => '<tr><td>' + esc(AG_LABEL[k] || k) + '</td><td class="mono">' + esc(cur[k] || '—') + '</td><td class="mono">' + esc(patch[k]) + '</td></tr>').join(''); }
/* 返回 {severity, fixes, retry, patch} 或 null */
async function agDiagnose(c, kind, text){
  if (!state.ai.has_key) { agNote(tr('没有配置 DeepSeek，无法自动诊断'), tr('请到“集群连接”页的文件浏览里查看日志，或先保存 API Key。')); return null; }
  const st = agStep(tr('读取最近的日志并请 DeepSeek 分析'));
  let r;
  try { r = await api('POST', '/api/agent/diagnose', {cluster:c.id, dir:genProjectDir(), kind, text:text || '', cfg:agCfgForAi()}); }
  catch (e) { st.fail(e.message); return null; }
  st.ok(r.logs.length ? tr('分析了这些日志：\n') + r.logs.join('\n') : tr('没有找到日志文件，只分析了上一步的输出。'));
  const patch = agMergePatch(r.fixes);
  agItem(r.severity === 'error' ? 'err' : '', '<div class="ag-h"><span class="pill ' + (r.severity === 'ok' ? 'ok' : r.severity === 'error' ? 'err' : 'warn') + '">' + ({ok:tr('正常'), warn:tr('注意'), error:tr('有问题')}[r.severity]) + tr('</span><b>诊断结果</b></div>') +
    '<p>' + esc(r.summary) + '</p>' + (r.cause ? tr('<p><b>可能原因：</b>') + esc(r.cause) + '</p>' : '') +
    (r.fixes.length ? '<ul>' + r.fixes.map(f => '<li>' + esc(f.text) + '</li>').join('') + '</ul>' : '') + tr('<p class="muted">以上是 AI 的判断，可能有误，请结合日志核对。</p>'));
  return {severity:r.severity, fixes:r.fixes, retry:r.retry, patch};
}
async function agDiagBtn(){
  if (!agBegin()) return;
  try {
    const c = agNeedCluster(); if (!c || !agNeedAi()) return;
    if (genErr() && !lastGen.cfg.project) { agErr(tr('先填写项目'), tr('表单里的项目名和工作目录是找日志的依据。')); return; }
    const d = await agDiagnose(c, 'pipeline', '');
    if (d && Object.keys(d.patch).length && await agAsk(tr('把建议的参数改动应用到表单？'), tr('<table class="ag-diff"><tr><th>字段</th><th>现在</th><th>改成</th></tr>') + agDiffRows(d.patch) + tr('</table><p>只改表单，不会重新提交。要重跑请用“一键跑完整流程”，或回到“流程生成”页手动部署。</p>'), tr('应用'))) {
      loadCfg(Object.assign(readCfg(), d.patch)); agNote(tr('已应用'), tr('表单已更新。'));
    }
  } finally { agEnd(); }
}

/* ---- ④ 端到端 ---- */
const agIds = out => { const s = new Set(); String(out).split('\n').forEach(l => { let m = l.match(/(?:作业号|Job ID):\s*(\d+)/) || l.match(/Submitted batch job (\d+)/) || l.match(/Your job (\d+)/) || l.match(/^(\d+)(?:\.\S+)?\s*$/); if (m) s.add(m[1]); }); return s; };
function agMyRows(out, ids){
  return String(out).split('\n').filter(l => { const t = l.trim(); if (!t) return false; const id = t.split(/[\s|]+/)[0].replace(/[._\[].*$/, ''); return ids.size ? ids.has(id) : /^\d+/.test(t); });
}
async function agDeploy(c){
  const dir = genProjectDir(), st = agStep(tr('上传并运行 deploy_kmeria.sh'));
  try {
    await api('POST', '/api/clusters/' + c.id + '/upload', {dir, files:[{name:'deploy_kmeria.sh', content:lastGen.deploy}]});
    const r = await api('POST', '/api/clusters/' + c.id + '/exec', {dir, action:'bash', script:'deploy_kmeria.sh'});
    if (r.exit_code === 0) { st.ok(fmtExec(r)); return {ok:true}; }
    st.fail(fmtExec(r)); return {ok:false, text:r.output};
  } catch (e) { st.fail(e.message); return {ok:false, text:e.message}; }
}
async function agSubmit(c){
  const dir = genProjectDir(), st = agStep(tr('提交作业链（submit_chain.sh --go）'));
  try {
    const r = await api('POST', '/api/clusters/' + c.id + '/exec', {dir, action:'bash', script:'submit_chain.sh', args:['--go']});
    if (r.exit_code === 0) { st.ok(fmtExec(r)); return {ok:true, ids:agIds(r.output)}; }
    st.fail(fmtExec(r)); return {ok:false, text:r.output};
  } catch (e) { st.fail(e.message); return {ok:false, text:e.message}; }
}
const AG_NOT_RESULT = /\.(sh|sbatch|job|pbs|qsub|err|out|log)$/;
async function agHasResult(c, dir){
  try { const r = await api('POST', '/api/clusters/' + c.id + '/ls', {path: dir}); return r.entries.filter(e => e.type === 'file' && !AG_NOT_RESULT.test(e.name) && e.size > 0).map(e => e.name); } catch (e) { return []; }
}
/* 监控到作业全部结束: 返回 {state:'done'|'failed'|'stopped', text} */
async function agMonitor(c, ids, resultDir, label){
  const st = agStep(tr('监控作业：') + label + tr('（每 ') + Math.max(1, Math.round(ag.pollMs / 1000)) + tr(' 秒查看一次）'));
  let gone = 0, n = 0;
  if (!ids.size) st.log(tr('没能从提交输出里读到作业号，改为监控你账号下的全部作业。'));
  while (!ag.stop) {
    let out = '';
    try { out = (await api('POST', '/api/clusters/' + c.id + '/jobs', {})).output || ''; } catch (e) { st.log(tr('查询失败：') + e.message + tr('，稍后重试')); await agSleep(ag.pollMs); continue; }
    const mine = agMyRows(out, ids); n++;
    st.log(tr('第 ') + n + tr(' 次查看，') + new Date().toLocaleTimeString() + '\n' + (mine.length ? mine.join('\n') : tr('队列里没有相关作业')));
    if (mine.some(l => /DependencyNeverSatisfied/.test(l))) { st.fail(); return {state:'failed', text:tr('有作业因为前一步失败而永远无法开始 (DependencyNeverSatisfied):\n') + mine.join('\n')}; }
    gone = mine.length ? 0 : gone + 1;
    if (gone >= 2) {
      const files = await agHasResult(c, resultDir);
      if (files.length) { st.ok(tr('作业都已结束，在 ') + resultDir + tr(' 里找到：') + files.slice(0, 6).join('、')); return {state:'done'}; }
      st.fail(); return {state:'failed', text:tr('作业都已结束，但 ') + resultDir + tr(' 里没有结果文件。')};
    }
    await agSleep(ag.pollMs);
  }
  st.skip(tr('已停止监控。作业仍在服务器上运行，可到“集群连接”页查看。')); return {state:'stopped'};
}
async function agPostGwas(c){
  pgSync();
  const errs = lastPg.checks.filter(x => x.level === 'err');
  if (errs.length) { agItem('err', tr('<div class="ag-h"><span class="pill err">还差 ') + errs.length + tr(' 项</span><b>Post-GWAS 设置未完成</b></div><ul>') + errs.slice(0, 5).map(x => '<li>' + esc(x.msg) + '</li>').join('') + tr('</ul><p><a href="#postgwas">到 Post-GWAS 页补充（通常是参考基因组 FASTA 路径）</a>，补完后可以直接在那一页提交。</p>')); return; }
  const o = lastPg.cfg;
  if (!await agAsk(tr('提交 Post-GWAS 定位与作图？'), tr('<p>阈值 p &lt; ') + esc(o.thresh) + tr('，BLAST 比对 + 画曼哈顿图，作为一个作业提交到集群。</p>'), tr('批准提交'))) return;
  const dir = pgOutDir(), st = agStep(tr('上传并提交 postgwas.sh')); let ids;
  try {
    await api('POST', '/api/clusters/' + c.id + '/upload', {dir, files:lastPg.files});
    const r = await api('POST', '/api/clusters/' + c.id + '/exec', {dir, action:'submit', script:'postgwas.sh'});
    if (r.exit_code !== 0) { st.fail(fmtExec(r)); await agDiagnose(c, 'postgwas', r.output); return; }
    st.ok(fmtExec(r)); ids = agIds(r.output);
  } catch (e) { st.fail(e.message); return; }
  const m = await agMonitor(c, ids, dir, 'Post-GWAS');
  if (m.state === 'failed') { await agDiagnose(c, 'postgwas', m.text); return; }
  if (m.state !== 'done') return;
  try {
    const t = await api('POST', '/api/clusters/' + c.id + '/tail', {path: dir + '/sigkmer.assoc.blastn.plot.txt', lines: 500});
    if (t.output.trim()) {
      $('pg_plotText').value = t.output; document.querySelector('input[name=pg_src][value=mine]').checked = true; updatePg();
      agNote(tr('结果已就绪'), tr('曼哈顿图（PDF/PNG）和 top_hits.tsv 在服务器的 ') + esc(dir) + tr('。定位结果的最后 500 行已填入 <a href="#postgwas">Post-GWAS 页</a>的预览，按你选的配色画出。'));
    }
  } catch (e) { agNote(tr('作业完成'), tr('结果在 ') + esc(dir) + tr('，读取预览数据失败：') + esc(e.message)); }
}
async function agAuto(){
  if (!agBegin()) return;
  try {
    const c = agNeedCluster(); if (!c) return;
    updateGen(); if (genErr()) { agCheckNote(); return; }
    agNote(tr('计划'), tr('① 检查环境（必要时安装，需批准）→ ② 部署（需批准）→ ③ 提交（需批准）→ ④ 监控 → 失败时诊断并提议修正（需批准，最多重试 2 次）→ ⑤ Post-GWAS（需批准）。'));
    if (!await agEnsureKmeria(c) || ag.stop) return;
    for (let attempt = 0; ; attempt++) {
      if (genErr()) { agCheckNote(); return; }
      if (!await agAsk(attempt ? tr('修正后重新部署') : tr('部署项目'), tr('<p>把 <code>deploy_kmeria.sh</code> 上传到 <code>') + esc(genProjectDir()) + tr('</code> 并运行：建目录、写输入文件、检查软件和数据、生成作业脚本。<b>不会提交作业。</b></p>'), tr('批准部署'))) return;
      let fail = null;
      const d = await agDeploy(c); if (!d.ok) fail = d.text;
      let ids = new Set();
      if (!fail) {
        if (!await agAsk(tr('提交作业'), tr('<p>运行 <code>submit_chain.sh --go</code>，按 ') + esc(lastGen.cfg.scheduler) + tr(' 顺序提交各步骤作业（后一步等前一步成功）。') + (attempt ? tr('<b>会重新提交整个流程，已有的中间结果是否复用取决于 KMERIA 本身。</b>') : '') + '</p>', tr('批准提交'))) return;
        const s = await agSubmit(c); if (!s.ok) fail = s.text; else ids = s.ids;
      }
      if (!fail) {
        const last = lastGen.cfg.step === 'all' ? 'asso' : lastGen.cfg.step;
        const m = await agMonitor(c, ids, joinPath(genProjectDir(), 'results/' + (STEP_DIR[last] || '05_association')), tr('主流程'));
        if (m.state === 'done') break; if (m.state === 'stopped') return; fail = m.text;
      }
      const dg = await agDiagnose(c, 'pipeline', fail);
      if (!dg) return;
      if (attempt >= 2) { agNote(tr('已自动重试 2 次仍未成功'), tr('请根据上面的诊断手动处理。')); return; }
      if (!Object.keys(dg.patch).length) { agNote(tr('没有可以自动应用的修改'), tr('需要你按建议处理（比如核对路径、环境），处理好后再点一次“一键跑完整流程”。')); return; }
      if (!await agAsk(tr('应用修正并重试（第 ') + (attempt + 1) + tr('/2 次）'), tr('<table class="ag-diff"><tr><th>字段</th><th>现在</th><th>改成</th></tr>') + agDiffRows(dg.patch) + '</table>', tr('应用并重试'))) return;
      loadCfg(Object.assign(readCfg(), dg.patch));
      if (genErr()) { agCheckNote(); return; }
    }
    agNote(tr('主流程已完成')); await agPostGwas(c);
  } finally { agEnd(); }
}

function renderAgSel(){
  const sel = $('ag_sel'); const keep = sel.value || state.activeId;
  sel.innerHTML = state.clusters.length ? state.clusters.map(c => '<option value="' + esc(c.id) + '">' + esc(c.name) + tr('（') + esc(c.username) + '@' + esc(c.host) + tr('）</option>')).join('') : tr('<option value="">（还没有连接）</option>');
  if (state.clusters.some(c => c.id === keep)) sel.value = keep;
}
function initAgent(){
  const pm = Number(LS.get('easykmera.pollMs')); if (pm >= 200) ag.pollMs = pm;
  routeHooks.assistant = () => { $('agStatic').hidden = state.backend; renderAgSel(); updateGen(); };
  $('ag_sel').addEventListener('change', () => { state.activeId = $('ag_sel').value; LS.set('easykmera.active', state.activeId); });
  $('ag_fill').addEventListener('click', agFill); $('ag_install').addEventListener('click', agInstall);
  $('ag_diag').addEventListener('click', agDiagBtn); $('ag_auto').addEventListener('click', agAuto);
  $('ag_stop').addEventListener('click', () => { ag.stop = true; if (ag.wake) ag.wake(); });
  $('ag_clear').addEventListener('click', () => { if (ag.busy) { toast(tr('任务进行中，不能清空')); return; } $('ag_log').innerHTML = tr('<p class="muted" id="ag_empty">还没有任务。左边选一个开始。</p>'); });
  $('ag_text').addEventListener('keydown', e => { if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); agFill(); } });
}

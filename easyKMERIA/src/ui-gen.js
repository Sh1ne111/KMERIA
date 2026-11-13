/* ===== 流程生成页 ===== */
const T = {
  deploy: trTpl($('tpl-deploy').textContent), submit: trTpl($('tpl-submit').textContent),
  postgwas: trTpl($('tpl-postgwas').textContent), plotR: trTpl($('tpl-plotr').textContent)
};
const GEN_FIELDS = Object.keys(DEFAULTS).filter(k => k !== 'dataMode');
let lastGen = null, genTimer = null;

function setVal(id, v){ const el = $(id); if (!el) return; if (el.type === 'checkbox') el.checked = !!v; else el.value = v; }
function getVal(id){ const el = $(id); if (!el) return ''; return el.type === 'checkbox' ? el.checked : el.value.trim(); }
function readCfg(){
  const cfg = {};
  GEN_FIELDS.forEach(k => { cfg[k] = (k === 'sampleText' || k === 'phenoText' || k === 'preCmd') ? $(k).value : getVal(k); });
  cfg.dataMode = document.querySelector('input[name=dataMode]:checked').value;
  return cfg;
}
function loadCfg(cfg){
  GEN_FIELDS.forEach(k => { if (k !== 'trait') setVal(k, cfg[k]); });
  document.querySelectorAll('input[name=dataMode]').forEach(r => { r.checked = r.value === cfg.dataMode; });
  refreshTraits(cfg.trait); applyMode(); updateGen();
}
function applyMode(){
  const mode = document.querySelector('input[name=dataMode]:checked').value;
  $('fastqDirBox').hidden = mode !== 'dir';
  if (mode === 'dir') {
    $('sampleLabel').textContent = tr('样本表');
    $('sampleText').placeholder = 'SS001\t31.2\nSS002\t28.7\nSS003\t35.4';
    $('sampleHint').innerHTML = tr('每行一个样本：<code>样本名</code> 空格或制表符 <code>测序深度</code>，可再加一列 <code>倍性</code>（混合倍性时才需要）。深度可以不写，会用下面的默认值。');
  } else {
    $('sampleLabel').textContent = tr('样本表（含数据路径）');
    $('sampleText').placeholder = 'SS001\t/data/a/SS001_1.fq.gz\t/data/a/SS001_2.fq.gz\t31.2\nSS002\t/data/b/SS002_1.fq.gz\t/data/b/SS002_2.fq.gz\t28.7';
    $('sampleHint').innerHTML = tr('每行一个样本：<code>样本名</code> <code>R1 路径</code> <code>R2 路径</code> <code>测序深度</code>，可再加一列 <code>倍性</code>。脚本会用软链接统一命名，原始数据不会被改动。');
  }
}
function refreshTraits(selected){
  const ph = parsePheno($('phenoText').value), sel = $('trait');
  const cur = selected != null ? String(selected) : sel.value;
  sel.innerHTML = '';
  if (!ph.traits.length) { const o = document.createElement('option'); o.value = '0'; o.textContent = tr('（还没有表型数据）'); sel.appendChild(o); return; }
  ph.traits.forEach((t, i) => { const o = document.createElement('option'); o.value = String(i); o.textContent = t + tr('（第 ') + (i + 1) + tr(' 列）'); sel.appendChild(o); });
  sel.value = (Number(cur) < ph.traits.length) ? cur : '0';
}
function checksHtml(res, title){
  const order = {err:0, warn:1, ok:2}, label = {err:tr('错误'), warn:tr('注意'), ok:tr('正常')};
  const items = res.checks.slice().sort((a, b) => order[a.level] - order[b.level]);
  const nErr = items.filter(i => i.level === 'err').length;
  return '<h3>' + (nErr ? tr('有 ') + nErr + tr(' 项错误需要先处理') : title) + '</h3>' +
    items.map(i => '<div class="chk"><span class="pill ' + i.level + '">' + label[i.level] + '</span><span>' + esc(i.msg) + '</span></div>').join('');
}
function renderGen(res){
  $('checks').innerHTML = checksHtml(res, tr('检查清单'));
  $('p-deploy').innerHTML = blk('deploy', 'deploy_kmeria.sh', res.deploy) +
    tr('<p class="note">脚本内部已经包含 <code>submit_chain.sh</code>，运行后会自动写出。整个脚本只做准备和预览，不会提交作业。</p>');
  $('p-inputs').innerHTML = blk('sample', 'sample.list', res.files.sampleList) + blk('depth', 'sample_depth.tsv', res.files.depthFile) + blk('pheno', 'phenotypes.txt', res.files.phenoFile) +
    tr('<p class="note">这三个文件已经包含在部署脚本里，这里单独列出是方便你核对内容。表型文件按样本列表的顺序排列，缺失的样本填 NA。</p>');
  $('p-cmd').innerHTML = blk('cmd', tr('kmeria_wrapper.pl 命令'), res.wrapperLit) +
    tr('<p class="note">这是部署脚本里调用官方 wrapper 的那一条命令，路径已经展开。想在 KMERIA 官方 Wiki 里对照参数时用它。</p>');
  $('p-install').innerHTML = blk('install', tr('安装命令'), res.install) +
    tr('<p class="note">来自官方 README 的预编译方式。需要 Linux、C/C++ 编译器和 GNU make，安装路径用的是左侧第 1 步填写的值。</p>');
  const s = res.summary;
  $('sampleSum').innerHTML = s.n ? tr('已识别 <b>') + s.n + tr('</b> 个样本') + (s.header ? tr('（第一行当作表头跳过）') : '') : '';
  const ph = parsePheno($('phenoText').value);
  $('phenoSum').innerHTML = ph.traits.length ? tr('已识别 <b>') + ph.rows.size + tr('</b> 行、<b>') + ph.traits.length + tr('</b> 个性状') + (ph.header ? tr('（含表头）') : tr('（没有表头，性状按 trait1、trait2 命名）')) : '';
}
function updateGen(){
  const cfg = readCfg();
  lastGen = generate(cfg, {deploy: T.deploy, submit: T.submit});
  lastGen.cfg = cfg;
  renderGen(lastGen);
  LS.set('easykmera.gen', cfg);
  const pc = $('p-cluster'); if (!pc.hidden) { if (pc._refresh) pc._refresh(); else renderGenRemote(); }
}
const scheduleGen = () => { clearTimeout(genTimer); genTimer = setTimeout(updateGen, 120); };
const genProjectDir = () => joinPath(lastGen.cfg.workDir || '/path/to/workdir', lastGen.cfg.project || 'kmeria_project');
const genErr = () => lastGen.checks.some(c => c.level === 'err') ? tr('检查清单里还有错误，请先处理再部署') : '';

function renderGenRemote(){
  remotePanel($('p-cluster'), {
    id:'rg',
    apply: c => {
      const cfg = readCfg();
      if (c.workdir) cfg.workDir = c.workdir; if (c.kmeria_home) cfg.kmeriaHome = c.kmeria_home;
      if (c.queue) cfg.queue = c.queue; cfg.scheduler = c.scheduler;
      loadCfg(cfg);
    },
    warn: c => c.scheduler !== lastGen.cfg.scheduler ? tr('注意：表单里的调度系统和这个连接不一致，提交作业前请先点“用它的设置填充表单”。') : '',
    logDir: () => genProjectDir() + '/results',
    steps: [
      {title:tr('1. 上传并运行部署脚本'), desc:tr('把 deploy_kmeria.sh 传到项目目录并运行：建目录、写输入文件、检查软件和数据、生成作业脚本。不会提交作业。'), label:tr('上传并运行'), primary:true, check: genErr,
       run: async c => {
         const dir = genProjectDir();
         await api('POST', '/api/clusters/' + c.id + '/upload', {dir, files:[{name:'deploy_kmeria.sh', content:lastGen.deploy}]});
         return fmtExec(await api('POST', '/api/clusters/' + c.id + '/exec', {dir, action:'bash', script:'deploy_kmeria.sh'}));
       }},
      {title:tr('2. 提交作业'), desc:tr('运行 submit_chain.sh --go，按顺序提交各步骤作业，后一步等前一步成功。请先完成第 1 步并核对预览。'), label:tr('提交作业'), arm:true, check: genErr,
       run: async c => fmtExec(await api('POST', '/api/clusters/' + c.id + '/exec', {dir:genProjectDir(), action:'bash', script:'submit_chain.sh', args:['--go']}))},
      {title:tr('3. 查看作业状态'), desc:tr('列出你账号下正在排队和运行的作业。'), label:tr('刷新状态'),
       run: async c => { const r = await api('POST', '/api/clusters/' + c.id + '/jobs', {}); return (r.output || '').trim() || tr('当前没有排队或运行中的作业。'); }}
    ]
  });
}

function initGenerate(){
  (function(){ const s = 'ATGCTTAGCCATGGCTAAGTCCGATTACGGA'; $('ribbon').innerHTML = s.split('').map(c => '<span class="nt-' + c + '">' + c + '</span>').join(''); })();
  $('form').addEventListener('input', e => {
    $('demo').hidden = true;
    if (e.target.id === 'phenoText') refreshTraits();
    if (e.target.name === 'dataMode') applyMode();
    scheduleGen();
  });
  $('form').addEventListener('change', e => { if (e.target.name === 'dataMode') applyMode(); scheduleGen(); });
  const bindFile = (inputId, targetId) => $(inputId).addEventListener('change', e => {
    const f = e.target.files && e.target.files[0]; if (!f) return;
    const r = new FileReader();
    r.onload = () => { $(targetId).value = String(r.result || ''); if (targetId === 'phenoText') refreshTraits(0); $('demo').hidden = true; updateGen(); };
    r.readAsText(f); e.target.value = '';
  });
  bindFile('fileSamples', 'sampleText'); bindFile('filePheno', 'phenoText');
  $('btnExample').addEventListener('click', () => { loadCfg(EXAMPLE); $('demo').hidden = false; });
  $('btnClear').addEventListener('click', () => { loadCfg(Object.assign({}, DEFAULTS)); $('demo').hidden = true; });
  bindTabs('tabs', 'p-', tab => { if (tab === 'cluster') renderGenRemote(); });
  const saved = LS.get('easykmera.gen');
  if (saved && saved.project) { loadCfg(Object.assign({}, DEFAULTS, saved)); $('demo').hidden = true; }
  else { loadCfg(EXAMPLE); $('demo').hidden = false; }
}

/* ===== Post-GWAS 页 ===== */
const PG_IDS = ['assocGlob','thresh','lineP','genomeFasta','blastDb','evalue','wordSize','blastThreads','excludeRe','useFai','genesOn','gff3','window','figW','figH','figTitle','pgMemory','pgTime','ptSize','hiSize','alpha','ymax','showLine','lineStyle','line2P','grid','stripChr','shape','dpi','format','labelTop','fontScale'];
let pgState = Object.assign({}, PG_DEFAULTS, {palette: JSON.parse(JSON.stringify(PALETTES[0]))});
let lastPg = null, pgTimer = null, demoCache = null;

function readPgForm(){
  const o = {};
  PG_IDS.forEach(k => { const el = $('pg_' + k); o[k] = el.type === 'checkbox' ? el.checked : el.value.trim(); });
  o.dbMode = document.querySelector('input[name=pg_dbMode]:checked').value;
  o.palette = pgState.palette;
  return o;
}
function writePgForm(o){
  PG_IDS.forEach(k => { const el = $('pg_' + k); if (el.type === 'checkbox') el.checked = o[k] === true || o[k] === 'true'; else el.value = o[k] == null ? '' : o[k]; });
  document.querySelectorAll('input[name=pg_dbMode]').forEach(r => { r.checked = r.value === o.dbMode; });
  pgState.palette = o.palette || pgState.palette;
}
function genContext(){
  const c = lastGen ? lastGen.cfg : readCfg();
  return {project:c.project, workDir:c.workDir, scheduler:c.scheduler, queue:c.queue, condaEnv:c.condaEnv, preCmd:c.preCmd, assoTool:c.assoTool};
}
function pgSync(){
  const g = genContext();
  if (!$('pg_assocGlob').value.trim() && g.workDir && g.project) $('pg_assocGlob').value = joinPath(joinPath(g.workDir, g.project), 'results/05_association/*.assoc.txt');
  updatePg();
}
const pgOutDir = () => lastPg.outDir;

function renderPresets(){
  const pk = o => JSON.stringify(['hi','line','bg','text','chr'].map(k => o[k]));
  const cur = pk(pgState.palette);
  $('pgPresets').innerHTML = PALETTES.map((p, i) => '<button type="button" class="preset" data-i="' + i + '" aria-pressed="' + (pk(p) === cur) + '"><span class="sw">' +
    p.chr.concat([p.hi]).map(c => '<i style="background:' + c + '"></i>').join('') + '</span>' + esc(p.name) + '</button>').join('');
  const p = pgState.palette;
  const items = p.chr.map((c, i) => ['chr' + i, tr('染色体色 ') + (i + 1), c]).concat([['hi',tr('高亮点'),p.hi],['line',tr('阈值线'),p.line],['bg',tr('背景'),p.bg],['text',tr('文字'),p.text]]);
  $('pgColors').innerHTML = items.map(([k, label, v]) => '<label><input type="color" data-k="' + k + '" value="' + esc(v) + '" aria-label="' + esc(label) + '">' + esc(label) + '</label>').join('');
}
function setPalette(p){ pgState.palette = p; renderPresets(); updatePg(); }

function updatePg(){
  const o = readPgForm();
  lastPg = generatePostGwas(o, genContext(), T);
  lastPg.cfg = o;
  $('pgChecks').innerHTML = checksHtml(lastPg, tr('检查清单'));
  $('pp-script').innerHTML = blk('pgscript', 'postgwas.sh', lastPg.script) + tr('<p class="note">会放在项目目录的 <code>postgwas</code> 子目录。头部已经带了集群资源设置，可以直接 <code>sbatch postgwas.sh</code>，也可以在登录节点用 <code>bash postgwas.sh</code> 运行。</p>');
  $('pp-rscript').innerHTML = blk('pgr', 'plot_manhattan_easy.R', lastPg.rscript) + tr('<p class="note">与 postgwas.sh 放在同一目录，只依赖 R 自带的基础绘图，不用装任何包。配色来自第 4 步的设置。</p>');
  $('pg_geneBox').hidden = !o.genesOn; $('pg_genomeBox').hidden = o.dbMode !== 'genome'; $('pg_dbBox').hidden = o.dbMode !== 'db';
  LS.set('easykmera.pg', Object.assign({}, o));
  drawPreview();
  const pc = $('pp-cluster'); if (!pc.hidden) { if (pc._refresh) pc._refresh(); else renderPgRemote(); }
}
const pgErr = () => lastPg.checks.some(c => c.level === 'err') ? tr('检查清单里还有错误，请先处理再部署') : '';
function renderPgRemote(){
  remotePanel($('pp-cluster'), {
    id:'rp',
    warn: c => { const g = genContext(); return c.scheduler !== g.scheduler ? tr('注意：脚本头部按“流程生成”里的调度系统（') + g.scheduler + tr('）生成，和这个连接（') + c.scheduler + tr('）不一致。') : ''; },
    logDir: () => pgOutDir(),
    steps: [
      {title:tr('1. 上传脚本'), desc:tr('把 postgwas.sh 和 plot_manhattan_easy.R 传到项目的 postgwas 目录。'), label:tr('上传'), check: pgErr,
       run: async c => { const r = await api('POST', '/api/clusters/' + c.id + '/upload', {dir:pgOutDir(), files:lastPg.files}); return tr('已上传：\n') + r.saved.join('\n'); }},
      {title:tr('2. 提交到集群'), desc:tr('用调度系统提交 postgwas.sh，BLAST 比对会在计算节点上运行。推荐。'), label:tr('提交作业'), primary:true, check: pgErr,
       run: async c => { await api('POST', '/api/clusters/' + c.id + '/upload', {dir:pgOutDir(), files:lastPg.files}); return fmtExec(await api('POST', '/api/clusters/' + c.id + '/exec', {dir:pgOutDir(), action:'submit', script:'postgwas.sh'})); }},
      {title:tr('或者：在登录节点直接运行'), desc:tr('只适合显著 k-mer 很少的情况，运行时间长会被管理员中断。'), label:tr('直接运行'), arm:true, check: pgErr,
       run: async c => { await api('POST', '/api/clusters/' + c.id + '/upload', {dir:pgOutDir(), files:lastPg.files}); return fmtExec(await api('POST', '/api/clusters/' + c.id + '/exec', {dir:pgOutDir(), action:'bash', script:'postgwas.sh'})); }},
      {title:tr('查看作业状态'), desc:tr('列出你账号下正在排队和运行的作业。'), label:tr('刷新状态'),
       run: async c => { const r = await api('POST', '/api/clusters/' + c.id + '/jobs', {}); return (r.output || '').trim() || tr('当前没有排队或运行中的作业。'); }}
    ]
  });
}

/* ---- 曼哈顿图预览 ---- */
function drawPreview(){
  const cv = $('pgCanvas'), W = cv.parentElement.clientWidth;
  const src = document.querySelector('input[name=pg_src]:checked').value;
  $('pg_mineBox').hidden = src !== 'mine';
  let rows, fai = '';
  if (src === 'demo') { demoCache = demoCache || demoPlotRows(); rows = demoCache.rows; fai = demoCache.fai; $('pg_previewNote').textContent = tr('这是随机生成的示例数据（虚线按 1e-6 画），只用来预览配色，不是你的结果。'); }
  else { rows = parsePlotTable($('pg_plotText').value); fai = $('pg_faiText').value; $('pg_previewNote').textContent = rows.length ? '' : tr('粘贴 sigkmer.assoc.blastn.plot.txt 的内容后，这里会画出你的结果。'); }
  const o = lastPg.cfg, lineP = src === 'demo' ? '1e-6' : (o.lineP || o.thresh), L = layoutManhattan(rows, fai, o.excludeRe, lineP);
  if (!W) return;
  drawManhattan(cv, L, o.palette, {w:Number(o.figW) || 12, h:Number(o.figH) || 4.5, title:o.figTitle, ptSize:Number(o.ptSize) || 0.7, hiSize:Number(o.hiSize) || 1, alpha:Number(o.alpha) || 1,
    ymax:Number(o.ymax) || 0, showLine:o.showLine, lineStyle:o.lineStyle, line2P:Number(o.line2P) || 0, grid:o.grid, stripChr:o.stripChr, shape:o.shape, labelTop:parseInt(o.labelTop, 10) || 0, fontScale:Number(o.fontScale) || 1});
  const nSig = L.pts.filter(p => p.sig).length;
  $('pgSum').textContent = L.pts.length ? L.pts.length + tr(' 个 k-mer，') + L.chrs.length + tr(' 条染色体，') + nSig + tr(' 个超过阈值线') : '';
  const top = L.pts.slice().sort((a, b) => a.p - b.p).slice(0, 8);
  $('pgTop').innerHTML = top.length ? tr('<tr><th>k-mer</th><th>染色体</th><th>位置</th><th>p 值</th></tr>') + top.map(p => '<tr><td>' + esc(p.id) + '</td><td>' + esc(p.chr) + '</td><td>' + p.pos.toLocaleString() + '</td><td>' + p.p.toExponential(2) + '</td></tr>').join('') : '';
}
function niceStep(max){ const raw = max / 5, pow = Math.pow(10, Math.floor(Math.log10(raw))); const f = raw / pow; return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10) * pow; }
function drawManhattan(cv, L, pal, o){
  const W = cv.parentElement.clientWidth, H = Math.round(W * o.h / o.w), dpr = Math.min(window.devicePixelRatio || 1, 2), fs = o.fontScale;
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.height = H + 'px';
  const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const font = "'IBM Plex Sans','Noto Sans SC',sans-serif";
  g.fillStyle = pal.bg; g.fillRect(0, 0, W, H);
  const m = {l:Math.round(40 + 10 * fs), r:14, t:o.title ? 20 + 12 * fs : 14, b:Math.round(30 + 12 * fs)}, pw = W - m.l - m.r, ph = H - m.t - m.b;
  if (!L.pts.length || !L.total) { g.fillStyle = pal.text; g.font = '14px ' + font; g.textAlign = 'center'; g.fillText(tr('还没有可画的数据'), W / 2, H / 2); return; }
  const thrY = -Math.log10(L.thr), l2 = o.line2P > 0 && o.line2P < 1 ? -Math.log10(o.line2P) : null;
  const ymax = o.ymax > 0 ? o.ymax : Math.max(Math.max.apply(null, L.pts.map(p => p.y)), thrY, l2 || 0) * 1.08;
  const X = x => m.l + x / L.total * pw, Y = y => m.t + ph - y / ymax * ph;
  g.strokeStyle = pal.text; g.fillStyle = pal.text; g.globalAlpha = 1; g.lineWidth = 1;
  g.font = (11 * fs) + 'px ' + font; g.textAlign = 'right'; g.textBaseline = 'middle';
  const step = niceStep(ymax);
  for (let v = 0; v <= ymax; v += step) { if (o.grid) { g.globalAlpha = .12; g.beginPath(); g.moveTo(m.l, Y(v)); g.lineTo(m.l + pw, Y(v)); g.stroke(); } g.globalAlpha = 1; g.fillText(String(Math.round(v * 100) / 100), m.l - 6, Y(v)); }
  const cols = pal.chr, vis = L.pts.filter(p => p.y <= ymax), r1 = 2.2 * o.ptSize / 0.7, r2 = 3.1 * o.hiSize;
  const mark = (x, y, r) => {
    g.beginPath();
    if (o.shape === 'triangle') { g.moveTo(x, y - r * 1.15); g.lineTo(x + r * 1.05, y + r * .8); g.lineTo(x - r * 1.05, y + r * .8); g.closePath(); }
    else if (o.shape === 'square') g.rect(x - r * .85, y - r * .85, r * 1.7, r * 1.7);
    else g.arc(x, y, r, 0, 6.2832);
    if (o.shape === 'hollow') { g.lineWidth = 1; g.stroke(); } else g.fill();
  };
  g.globalAlpha = o.alpha;
  vis.filter(p => !p.sig).forEach(p => { g.fillStyle = g.strokeStyle = cols[p.ci % cols.length]; mark(X(p.x), Y(p.y), r1); });
  vis.filter(p => p.sig).forEach(p => { g.fillStyle = g.strokeStyle = pal.hi; mark(X(p.x), Y(p.y), r2); });
  g.globalAlpha = 1;
  g.strokeStyle = pal.line;
  if (o.showLine) { g.setLineDash(o.lineStyle === 'solid' ? [] : o.lineStyle === 'dotted' ? [2, 3] : [6, 4]); g.beginPath(); g.moveTo(m.l, Y(thrY)); g.lineTo(m.l + pw, Y(thrY)); g.stroke(); }
  if (l2 != null && l2 <= ymax) { g.setLineDash([2, 3]); g.beginPath(); g.moveTo(m.l, Y(l2)); g.lineTo(m.l + pw, Y(l2)); g.stroke(); }
  g.setLineDash([]);
  g.strokeStyle = pal.text; g.globalAlpha = .7; g.strokeRect(m.l, m.t, pw, ph); g.globalAlpha = 1;
  g.fillStyle = pal.text; g.textAlign = 'center'; g.textBaseline = 'top'; g.font = (11 * fs) + 'px ' + font;
  const every = pw / L.chrs.length < 24 * fs ? 2 : 1;
  L.chrs.forEach((c, i) => { if (i % every === 0) g.fillText(o.stripChr ? String(c.chr).replace(/^[Cc]hr/, '') : String(c.chr), X(c.off + c.len / 2), m.t + ph + 7); });
  g.font = (12 * fs) + 'px ' + font; g.fillText('Chromosome', m.l + pw / 2, m.t + ph + 8 + 14 * fs);
  g.save(); g.translate(13, m.t + ph / 2); g.rotate(-Math.PI / 2); g.textBaseline = 'top'; g.fillText('−log10(P)', 0, 0); g.restore();
  if (o.labelTop > 0) {
    g.font = (9 * fs) + 'px ' + font; g.textBaseline = 'bottom'; g.textAlign = 'center';
    vis.slice().sort((a, b) => a.p - b.p).slice(0, o.labelTop).forEach(p => g.fillText(p.id, X(p.x), Y(p.y) - 6));
  }
  if (o.title) { g.font = '600 ' + (14 * fs) + 'px ' + font; g.textBaseline = 'top'; g.fillText(o.title, W / 2, 9); }
}

const PG_AI_LABEL = {figTitle:tr('图标题'), figW:tr('图宽'), figH:tr('图高'), lineP:tr('阈值线 p 值'), ptSize:tr('普通点大小'), hiSize:tr('高亮点大小'), alpha:tr('不透明度'), ymax:tr('纵轴上限'), showLine:tr('画阈值线'),
  lineStyle:tr('阈值线线型'), line2P:tr('第二条阈值线'), grid:tr('横向网格'), stripChr:tr('去掉 Chr 前缀'), shape:tr('点形状'), dpi:tr('分辨率'), format:tr('输出格式'), labelTop:tr('标注前 N 个'), fontScale:tr('文字缩放'), palette:tr('配色')};
let pgUndo = null;
const pgStyleKeys = () => Object.keys(PG_AI_LABEL).filter(k => k !== 'palette');
async function pgAiGo(){
  const msg = $('pg_aiMsg'), want = $('pg_aiPrompt').value.trim();
  if (!state.backend) { msg.textContent = tr('需要运行 easyKMERIA 程序并在“AI 设置”里保存 DeepSeek API Key。'); return; }
  if (!want) { msg.textContent = tr('先写一句想怎么改。'); return; }
  $('pg_aiGo').disabled = true; msg.textContent = tr('正在请求 DeepSeek…');
  try {
    const cur = readPgForm(), nChr = demoCache && document.querySelector('input[name=pg_src]:checked').value === 'demo' ? 20 : (lastPg && lastPg.chrCount) || 0;
    const send = {palette:cur.palette}; pgStyleKeys().forEach(k => { send[k] = cur[k]; });
    const r = await api('POST', '/api/ai/plot', {prompt:want, current:send, n_chr:nChr});
    const before = JSON.parse(JSON.stringify(cur)), after = Object.assign({}, cur), changed = [];
    Object.keys(r.options).forEach(k => { if (String(after[k]) !== String(r.options[k])) { changed.push(tr(PG_AI_LABEL[k] || k) + ' ' + (before[k] === '' ? tr('自动') : String(before[k])) + ' → ' + (r.options[k] === '' ? tr('自动') : String(r.options[k]))); after[k] = r.options[k]; } });
    if (Object.keys(r.palette).length) { const p = Object.assign({}, before.palette, r.palette); if (r.palette.chr) p.chr = r.palette.chr; after.palette = p; changed.push(tr('配色') + (r.palette.name ? tr('：') + r.palette.name : '')); }
    let html = '';
    if (changed.length) { pgUndo = before; writePgForm(after); $('pg_aiUndo').hidden = false; renderPresets(); updatePg(); html += '<b>' + esc(tr('已修改')) + tr('</b>：') + changed.map(esc).join(tr('；')) + '。'; }
    else html += esc(tr('没有需要修改的参数。'));
    if (r.notes) html += ' ' + esc(r.notes);
    if (r.unsupported.length) html += '<br><b>' + esc(tr('这些现在做不到')) + tr('</b>：') + r.unsupported.map(esc).join(tr('；')) + '。';
    msg.innerHTML = html;
  } catch (e) { msg.textContent = e.message + (/API Key/.test(e.message) ? tr('（到“AI 设置”里保存）') : ''); }
  finally { $('pg_aiGo').disabled = false; }
}
function pgAiUndo(){
  if (!pgUndo) return;
  writePgForm(pgUndo); pgUndo = null; $('pg_aiUndo').hidden = true; $('pg_aiMsg').textContent = tr('已撤销。'); renderPresets(); updatePg();
}

function initPostGwas(){
  const saved = LS.get('easykmera.pg');
  if (saved) pgState = Object.assign({}, PG_DEFAULTS, saved, {palette: saved.palette || pgState.palette});
  writePgForm(pgState); renderPresets();
  $('pgForm').addEventListener('input', e => { if (e.target.dataset && e.target.dataset.k) return; clearTimeout(pgTimer); pgTimer = setTimeout(updatePg, 120); });
  $('pgForm').addEventListener('change', e => { if (e.target.name === 'pg_dbMode') updatePg(); });
  $('pgPresets').addEventListener('click', e => { const b = e.target.closest('.preset'); if (b) setPalette(JSON.parse(JSON.stringify(PALETTES[Number(b.dataset.i)]))); });
  $('pgColors').addEventListener('input', e => {
    const k = e.target.dataset.k; if (!k) return; const p = JSON.parse(JSON.stringify(pgState.palette));
    if (k.indexOf('chr') === 0) p.chr[Number(k.slice(3))] = e.target.value; else p[k] = e.target.value;
    p.name = tr('自定义'); pgState.palette = p; updatePg(); document.querySelectorAll('.preset').forEach(b => b.setAttribute('aria-pressed', 'false'));
  });
  $('pg_aiGo').addEventListener('click', pgAiGo); $('pg_aiUndo').addEventListener('click', pgAiUndo);
  $('pg_aiPrompt').addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); pgAiGo(); } });
  document.querySelectorAll('input[name=pg_src]').forEach(r => r.addEventListener('change', drawPreview));
  $('pg_plotText').addEventListener('input', drawPreview); $('pg_faiText').addEventListener('input', drawPreview);
  bindTabs('pgTabs', 'pp-', tab => { if (tab === 'cluster') renderPgRemote(); if (tab === 'preview') drawPreview(); });
  if (window.ResizeObserver) new ResizeObserver(() => { if (lastPg && !$('v-postgwas').hidden && !$('pp-preview').hidden) drawPreview(); }).observe($('pgCanvas').parentElement);
  $('pgSavePng').addEventListener('click', () => { $('pgCanvas').toBlob(b => { const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = 'manhattan_preview.png'; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }); });
  routeHooks.postgwas = pgSync;
  if (!$('pg_aiGo').dataset.init) { $('pg_aiGo').dataset.init = '1'; }
}

/* ===== 启动 ===== */
(async function boot(){
  i18nStatic(); initChrome();
  initGenerate(); initClusters(); initAi(); initAgent(); initPostGwas();
  $('btnCite').addEventListener('click', () => copyText($('citeText').textContent.trim(), $('btnCite'), tr('复制引用格式')));
  state.backend = await detectBackend();
  state.activeId = LS.get('easykmera.active') || '';
  const chip = $('modeChip');
  chip.textContent = state.backend ? tr('已连接程序') : tr('只读预览'); chip.classList.toggle('on', state.backend);
  chip.title = state.backend ? tr('已连接本机的 easyKMERIA 程序，可以登录服务器和使用 AI 配色') : tr('没有连接到 easyKMERIA 程序：脚本生成和预览可用，登录服务器与 AI 配色需要在本机运行程序');
  $('staticNotice').hidden = state.backend; $('pgSavePng').hidden = !state.backend;
  if (state.backend) { try { state.clusters = (await api('GET', '/api/clusters')).clusters; } catch (e) {} try { Object.assign(state.ai, await api('GET', '/api/ai/settings')); } catch (e) {} }
  route();
})();

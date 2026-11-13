const $ = id => document.getElementById(id);
const esc = s => String(s == null ? '' : s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
const PLACEHOLDER = '__EASYKMERA_' + 'TOKEN__';
const state = {backend:false, token:'', dataDir:'', clusters:[], activeId:'', ai:{has_key:false, base_url:'', model:''}};
const LS = {
  get(k){ try { return JSON.parse(localStorage.getItem(k) || 'null'); } catch (e) { return null; } },
  set(k, v){ try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) {} }
};
const ROUTES = ['home','generate','clusters','postgwas','assistant','ai'];
const routeHooks = {};

let toastTimer = null;
function toast(msg){
  const el = $('toast'); el.textContent = msg; el.hidden = false;
  clearTimeout(toastTimer); toastTimer = setTimeout(() => { el.hidden = true; }, 2600);
}
function copyText(text, btn, idle){
  const done = ok => { if (!btn) { toast(ok ? tr('已复制') : tr('复制失败，请手动选中后按 Ctrl+C')); return; }
    const old = idle || btn.textContent; btn.textContent = ok ? tr('已复制') : tr('请按 Ctrl+C'); setTimeout(() => { btn.textContent = old; }, 1600); };
  const fallback = () => {
    const ta = document.createElement('textarea'); ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
    document.body.appendChild(ta); ta.select(); let ok = false; try { ok = document.execCommand('copy'); } catch (e) {}
    ta.remove(); done(ok);
  };
  if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(() => done(true), fallback); else fallback();
}

/* 代码块 */
const blocks = {};
function hl(text){
  return esc(text).split('\n').map(line => /^\s*#/.test(line) ? '<span class="c">' + line + '</span>'
    : line.replace(/(^|\s)(--?[A-Za-z][\w-]*)/g, '$1<span class="f">$2</span>')).join('\n');
}
function blk(id, title, text){
  blocks[id] = text;
  return '<div class="blk"><div class="blk-h"><span class="fn">' + esc(title) + '</span><button type="button" class="copy" data-id="' + id + tr('">复制</button></div>') +
    '<pre class="code" tabindex="0"><code>' + hl(text) + '</code></pre></div>';
}
document.addEventListener('click', e => {
  const b = e.target.closest('.copy'); if (!b) return;
  copyText(blocks[b.dataset.id] || '', b, tr('复制'));
});

/* 标签页 */
function bindTabs(tabsId, prefix, onShow){
  const el = $(tabsId);
  el.addEventListener('click', e => {
    const b = e.target.closest('.tab'); if (!b) return;
    el.querySelectorAll('.tab').forEach(x => { x.setAttribute('aria-selected', String(x === b)); $(prefix + x.dataset.tab).hidden = x !== b; });
    if (onShow) onShow(b.dataset.tab);
  });
}

/* 弹窗 */
function openModal(id){ const m = $(id); m.hidden = false; const f = m.querySelector('input:not([type=hidden]):not([hidden])'); if (f) setTimeout(() => f.focus(), 30); }
function closeModal(id){ $(id).hidden = true; }
document.addEventListener('click', e => {
  if (e.target.closest('[data-close]')) { const m = e.target.closest('.modal'); if (m) m.hidden = true; }
  else if (e.target.classList && e.target.classList.contains('modal')) e.target.hidden = true;
});
document.addEventListener('keydown', e => { if (e.key === 'Escape') document.querySelectorAll('.modal').forEach(m => { m.hidden = true; }); });

/* 后端接口 */
async function api(method, path, body){
  let r;
  try {
    r = await fetch(path, {method, headers:{'Content-Type':'application/json', 'X-Token':state.token, 'X-Lang':LANG}, body: body === undefined ? undefined : JSON.stringify(body)});
  } catch (e) { throw new Error(tr('连不上 easyKMERIA 程序，请确认它还在运行')); }
  let j = null; try { j = await r.json(); } catch (e) {}
  if (!r.ok || !j || j.ok === false) throw new Error((j && j.error ? tr(j.error) : '') || (tr('请求失败（') + r.status + tr('）')));
  return j;
}
async function detectBackend(){
  const m = document.querySelector('meta[name="easykmera-token"]');
  const tok = m ? m.getAttribute('content') : '';
  if (!tok || tok === PLACEHOLDER) return false;
  try {
    const r = await fetch('/api/health'); const j = await r.json();
    if (!j || !j.ok) return false;
    state.token = tok; state.dataDir = j.data_dir || ''; state.version = j.version; return true;
  } catch (e) { return false; }
}

/* 路由 */
function route(){
  let h = (location.hash || '#home').slice(1);
  if (ROUTES.indexOf(h) < 0) h = 'home';
  ROUTES.forEach(r => { $('v-' + r).hidden = r !== h; });
  document.querySelectorAll('#navTabs a').forEach(a => { if (a.dataset.route === h) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current'); });
  if (routeHooks[h]) routeHooks[h]();
  window.scrollTo(0, 0);
}
window.addEventListener('hashchange', route);

/* 远程操作面板: 选连接, 上传, 运行, 看作业和日志 */
function remotePanel(el, spec){
  el._refresh = null;
  if (!state.backend) {
    el.innerHTML = tr('<div class="callout">这里需要运行 easyKMERIA 程序才能连接服务器。在你的电脑上运行 <code>python app.py</code> 后，用浏览器打开 <code>http://127.0.0.1:8765</code> 即可。</div>'); return;
  }
  if (!state.clusters.length) {
    el.innerHTML = tr('<p class="note">还没有保存的服务器连接。</p><p><a class="btn primary" href="#clusters">去添加连接</a></p>'); return;
  }
  if (!state.clusters.some(c => c.id === state.activeId)) state.activeId = state.clusters[0].id;
  const id = spec.id;
  const opts = state.clusters.map(c => '<option value="' + esc(c.id) + '"' + (c.id === state.activeId ? ' selected' : '') + '>' + esc(c.name) + tr('（') + esc(c.username) + '@' + esc(c.host) + tr('）</option>')).join('');
  el.innerHTML = '<div class="remote">' +
    '<div class="bar"><label class="lab" for="' + id + tr('_sel">服务器</label><select id="') + id + '_sel">' + opts + '</select>' +
    (spec.apply ? tr('<button type="button" class="btn" data-act="apply">用它的设置填充表单</button>') : '') + '</div>' +
    '<div class="note" id="' + id + '_warn"></div>' +
    '<div class="rsteps">' + spec.steps.map((s, i) => '<div class="rstep"><div class="t"><b>' + esc(s.title) + '</b><span class="hint">' + esc(s.desc) + '</span></div>' +
      '<button type="button" class="btn' + (s.primary ? ' primary' : '') + '" data-step="' + i + '">' + esc(s.label) + '</button></div>').join('') + '</div>' +
    '<pre class="rlog" id="' + id + '_log" hidden></pre>' +
    '<div class="field"><label class="lab" for="' + id + tr('_dir">浏览服务器上的目录和日志</label>') +
    '<div class="ai-row"><input type="text" id="' + id + '_dir" spellcheck="false" value="' + esc(spec.logDir()) + tr('"><button type="button" class="btn" data-act="ls">浏览</button></div>') +
    '<div class="filelist" id="' + id + '_files"></div></div></div>';
  const cur = () => state.clusters.find(c => c.id === $(id + '_sel').value);
  const log = $(id + '_log');
  const show = txt => { log.hidden = false; log.textContent = txt; log.scrollTop = log.scrollHeight; };
  const warn = () => { const c = cur(); $(id + '_warn').textContent = (spec.warn && c) ? spec.warn(c) : ''; };
  warn(); el._refresh = warn;
  $(id + '_sel').addEventListener('change', () => { state.activeId = $(id + '_sel').value; LS.set('easykmera.active', state.activeId); warn(); });
  async function ls(path){
    const c = cur(); if (!c) return;
    try {
      const r = await api('POST', '/api/clusters/' + c.id + '/ls', {path});
      $(id + '_dir').value = r.path;
      const up = r.path.replace(/\/[^\/]+\/?$/, '') || '/';
      $(id + '_files').innerHTML = '<button type="button" data-ls="' + esc(up) + '">../</button>' + r.entries.map(en =>
        '<button type="button" data-' + (en.type === 'dir' ? 'ls' : 'tail') + '="' + esc(r.path.replace(/\/$/, '') + '/' + en.name) + '">' + esc(en.name) + (en.type === 'dir' ? '/' : '') + '</button>').join('');
    } catch (e) { show(tr('浏览失败：') + e.message); }
  }
  el.onclick = async e => {
    const lsb = e.target.closest('[data-ls]'), tb = e.target.closest('[data-tail]'), step = e.target.closest('[data-step]'), act = e.target.closest('[data-act]');
    const c = cur();
    if (lsb) return ls(lsb.dataset.ls);
    if (tb) { try { const r = await api('POST', '/api/clusters/' + c.id + '/tail', {path: tb.dataset.tail, lines: 150}); show('— ' + tb.dataset.tail + tr(' 的最后 150 行 —\n') + r.output); } catch (er) { show(tr('读取失败：') + er.message); } return; }
    if (act && act.dataset.act === 'ls') return ls($(id + '_dir').value.trim());
    if (act && act.dataset.act === 'apply') { spec.apply(c); toast(tr('已填入这个连接的设置')); warn(); return; }
    if (step) {
      const s = spec.steps[Number(step.dataset.step)];
      if (s.arm && step.dataset.armed !== '1') {
        step.dataset.armed = '1'; const old = step.textContent; step.textContent = tr('再点一次确认'); step.classList.add('armed');
        setTimeout(() => { step.dataset.armed = ''; step.textContent = old; step.classList.remove('armed'); }, 5000); return;
      }
      if (s.arm) { step.dataset.armed = ''; step.classList.remove('armed'); step.textContent = s.label; }
      const pre = s.check ? s.check() : ''; if (pre) { toast(pre); return; }
      step.disabled = true; show(tr('运行中，请稍候…'));
      try { show(await s.run(c)); } catch (er) { show(tr('失败：') + er.message); } finally { step.disabled = false; }
    }
  };
}
const fmtExec = r => tr('退出码 ') + r.exit_code + (r.exit_code === 0 ? tr('（成功）') : tr('（有错误）')) + '\n\n' + (r.output || '').trim();

/* ===== 外壳: 语言和主题 ===== */
const THEMES = ['auto', 'light', 'dark'];
function applyTheme(m){
  const root = document.documentElement;
  if (m === 'auto') root.removeAttribute('data-theme'); else root.setAttribute('data-theme', m);
  $('themeLabel').textContent = tr({auto:tr('自动'), light:tr('浅色'), dark:tr('深色')}[m]);
}
function initChrome(){
  document.querySelectorAll('.seg [data-lang]').forEach(b => {
    b.setAttribute('aria-pressed', String(b.dataset.lang === LANG));
    b.addEventListener('click', () => {
      if (b.dataset.lang === LANG) return;
      if (typeof ag !== 'undefined' && ag.busy && !confirm(tr('AI 助手还有任务在进行，切换语言会重新加载页面并中断监控。确定切换吗？'))) return;
      setLang(b.dataset.lang);
    });
  });
  let mode = LS.get('easykmera.theme'); if (THEMES.indexOf(mode) < 0) mode = 'auto';
  applyTheme(mode);
  $('themeBtn').addEventListener('click', () => { mode = THEMES[(THEMES.indexOf(mode) + 1) % THEMES.length]; LS.set('easykmera.theme', mode); applyTheme(mode); });
}

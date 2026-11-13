/* ===== 中英文界面 =====
   语言在页面加载时确定 (读 localStorage, 没有就看浏览器语言); 切换语言后重新加载页面。
   代码里所有面向用户的文字都写成 tr('中文原文'); 英文翻译放在 i18n/en.json, 由 build.py 嵌入。 */
const LANG = (() => {
  let v = null; try { v = localStorage.getItem('easykmera.lang'); } catch (e) {}
  if (v === 'zh' || v === 'en') return v;
  return /^zh/i.test((navigator.languages && navigator.languages[0]) || navigator.language || '') ? 'zh' : 'en';
})();
const I18N_EN = /*I18N_EN*/{}/*END*/;
const I18N_PAT = [];
Object.keys(I18N_EN).forEach(k => {
  if (k.indexOf('{*}') < 0) return;
  I18N_PAT.push([new RegExp('^' + k.split('{*}').map(x => x.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')).join('([\\s\\S]*?)') + '$'), I18N_EN[k]]);
});
const CJK = /[㐀-鿿＀-￯]/;
function tr(s){
  if (LANG === 'zh' || typeof s !== 'string' || !CJK.test(s)) return s;
  const v = I18N_EN[s]; if (v !== undefined) return v;
  for (let i = 0; i < I18N_PAT.length; i++) {
    const m = s.match(I18N_PAT[i][0]);
    if (m) { let n = 1; return I18N_PAT[i][1].replace(/\{\*\}/g, () => m[n++] || ''); }
  }
  return s;
}
/* 脚本模板: 按整行翻译 (代码部分原样保留在译文里) */
function trTpl(text){
  if (LANG === 'zh') return text;
  return text.split('\n').map(l => (CJK.test(l) && I18N_EN[l.replace(/\s+$/, '')] !== undefined) ? I18N_EN[l.replace(/\s+$/, '')] : l).join('\n');
}
/* 页面里写死的 HTML: 逐个文本节点 / 属性 / 含标签的整段翻译 */
const I18N_ATTRS = ['placeholder', 'title', 'aria-label', 'alt'];
const i18nKey = s => s.replace(/\s+/g, ' ').trim();
function i18nWalk(el, visit){
  const tag = el.tagName;
  if (tag === 'SCRIPT' || tag === 'STYLE') return;
  I18N_ATTRS.forEach(a => { const v = el.getAttribute && el.getAttribute(a); if (v && CJK.test(v)) visit('attr', el, a, v); });
  if (tag === 'PRE' || tag === 'TEXTAREA') return;
  const kids = Array.prototype.slice.call(el.childNodes);
  const mixed = el.children.length > 0 && kids.some(n => n.nodeType === 3 && CJK.test(n.nodeValue));
  if (mixed && visit('html', el, null, i18nKey(el.innerHTML))) return;
  kids.forEach(n => {
    if (n.nodeType === 3) { if (CJK.test(n.nodeValue)) visit('text', n, null, i18nKey(n.nodeValue)); }
    else if (n.nodeType === 1) i18nWalk(n, visit);
  });
}
function i18nStatic(){
  document.documentElement.lang = LANG === 'en' ? 'en' : 'zh-CN';
  if (LANG !== 'en') return;
  i18nWalk(document.querySelector('.app'), (kind, node, attr, key) => {
    const v = tr(key); if (v === key) return false;
    if (kind === 'attr') node.setAttribute(attr, v);
    else if (kind === 'html') node.innerHTML = v;
    else { const m = node.nodeValue.match(/^(\s*)[\s\S]*?(\s*)$/); node.nodeValue = m[1] + v + m[2]; }
    return true;
  });
}
function setLang(l){
  if (l === LANG) return;
  try { localStorage.setItem('easykmera.lang', l); } catch (e) {}
  location.reload();
}

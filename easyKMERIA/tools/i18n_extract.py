#!/usr/bin/env python3
"""收集所有需要翻译的中文: 页面静态 HTML / JS 字符串 / 脚本模板行 / 后端提示。输出 src/i18n/todo.json
用法: python tools/i18n_extract.py"""
import ast, json, re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
S = ROOT / 'src'
CJK = re.compile(r'[㐀-鿿＀-￯]')
items = {}   # key -> {"kind":..., "ctx":...}
def add(key, kind, ctx=''):
    if key and CJK.search(key) and key not in items: items[key] = {'kind': kind, 'ctx': ctx}

# 1. 页面静态 HTML (用浏览器里同一套遍历逻辑, 保证键和运行时一致)
from playwright.sync_api import sync_playwright
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(locale='zh-CN')
    pg.goto((ROOT / 'static' / 'index.html').as_uri()); pg.wait_for_timeout(600)
    units = pg.evaluate("""() => { const out = [];
      i18nWalk(document.querySelector('.app'), (kind, node, attr, key) => { out.push([kind, key]); return kind === 'html'; });
      return out; }""")
    b.close()
for kind, key in units: add(key, 'html' if kind == 'html' else 'html-text' if kind == 'text' else 'html-attr')

# 2. JS 字符串
for x in json.loads((S / 'i18n' / 'js_literals.json').read_text(encoding='utf8')):
    f, ln = x['where'][0].split(':'); line = (S / f).read_text(encoding='utf8').split('\n')[int(ln) - 1].strip()
    add(x['zh'], 'js', f'{f}:{ln}  ' + line[:260])

# 3. 脚本模板: 整行
for f in ['tpl-deploy.txt', 'tpl-submit.txt', 'tpl-postgwas.txt', 'plot.R.txt', 'env-activate.sh']:
    for line in (S / f).read_text(encoding='utf8').split('\n'):
        if CJK.search(line): add(line.rstrip(), 'script-line', f)

# 4. 后端提示 (f-string 里的变量写成 {*})
src = (ROOT / 'app.py').read_text(encoding='utf8'); tree = ast.parse(src); lines = src.split('\n')
SKIP = {'PALETTE_SYSTEM', 'FILL_SYSTEM', 'DIAG_SYSTEM', 'PLOT_SYSTEM', 'FORM_HELP', 'PLOT_HELP'}
skip_nodes = set()
for node in ast.walk(tree):
    if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in SKIP for t in node.targets):
        for n in ast.walk(node): skip_nodes.add(id(n))
    if isinstance(node, (ast.FunctionDef, ast.ClassDef, ast.Module)) and node.body and isinstance(node.body[0], ast.Expr) and isinstance(getattr(node.body[0], 'value', None), ast.Constant):
        skip_nodes.add(id(node.body[0].value))
for node in ast.walk(tree):
    if id(node) in skip_nodes: continue
    if isinstance(node, ast.JoinedStr):
        parts = [v.value if isinstance(v, ast.Constant) else '{*}' for v in node.values]
        add(''.join(parts), 'backend', lines[node.lineno - 1].strip()[:200])
        for v in node.values: skip_nodes.add(id(v))
    elif isinstance(node, ast.Constant) and isinstance(node.value, str):
        add(node.value, 'backend', lines[node.lineno - 1].strip()[:200])
(S / 'i18n' / 'todo.json').write_text(json.dumps(items, ensure_ascii=False, indent=1), encoding='utf8')
from collections import Counter
print(len(items), dict(Counter(v['kind'] for v in items.values())))

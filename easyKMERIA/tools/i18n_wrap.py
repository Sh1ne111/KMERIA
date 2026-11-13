#!/usr/bin/env python3
"""把 JS 源码里含中文的字符串字面量包成 tr('...'), 并导出这些字符串。只需运行一次; 之后新增文字直接写 tr('中文')。
用法: python tools/i18n_wrap.py [--dry]   输出 src/i18n/js_literals.json"""
import json, re, sys
from pathlib import Path
S = Path(__file__).resolve().parent.parent / 'src'
FILES = ['core.js', 'postgwas.js', 'ui-common.js', 'ui-gen.js', 'ui-clusters.js', 'ui-agent.js', 'ui-pg.js']
CJK = re.compile(r'[㐀-鿿＀-￯]')
REGEX_PREV = set('(,=:[!&|?{};+-*%<>~^') 

def unescape(body):
    out, i = [], 0
    while i < len(body):
        c = body[i]
        if c == '\\' and i + 1 < len(body):
            n = body[i + 1]
            if n == 'u' and re.match(r'[0-9a-fA-F]{4}', body[i + 2:i + 6]): out.append(chr(int(body[i + 2:i + 6], 16))); i += 6; continue
            if n == 'x' and re.match(r'[0-9a-fA-F]{2}', body[i + 2:i + 4]): out.append(chr(int(body[i + 2:i + 4], 16))); i += 4; continue
            out.append({'n': '\n', 't': '\t', 'r': '\r', '0': '\0'}.get(n, n)); i += 2; continue
        out.append(c); i += 1
    return ''.join(out)

def lex(src):
    """返回 [(kind, text)]; kind: str / other"""
    toks, i, n = [], 0, len(src)
    last = ''   # 上一个有意义的字符, 用来判断 / 是除号还是正则
    buf = []
    def flush():
        if buf: toks.append(('o', ''.join(buf))); buf.clear()
    while i < n:
        c = src[i]
        if c == '/' and src[i:i + 2] == '//':
            j = src.find('\n', i); j = n if j < 0 else j; buf.append(src[i:j]); i = j; continue
        if c == '/' and src[i:i + 2] == '/*':
            j = src.find('*/', i + 2); j = n if j < 0 else j + 2; buf.append(src[i:j]); i = j; continue
        if c in '\'"':
            j = i + 1
            while j < n and src[j] != c:
                j += 2 if src[j] == '\\' else 1
            flush(); toks.append(('s', src[i:j + 1])); i = j + 1; last = c; continue
        if c == '`':
            j = i + 1
            while j < n and src[j] != '`': j += 2 if src[j] == '\\' else 1
            buf.append(src[i:j + 1]); i = j + 1; last = '`'; continue
        if c == '/':
            if last == '' or last in REGEX_PREV or re.search(r'(return|typeof|case)$', ''.join(buf)[-8:]):
                j, incls = i + 1, False
                while j < n:
                    if src[j] == '\\': j += 2; continue
                    if src[j] == '[': incls = True
                    elif src[j] == ']': incls = False
                    elif src[j] == '/' and not incls: break
                    j += 1
                j += 1
                while j < n and src[j].isalpha(): j += 1
                buf.append(src[i:j]); i = j; last = ')'; continue
        buf.append(c)
        if not c.isspace(): last = c
        i += 1
    flush()
    return toks

def main(dry):
    lits = []
    for f in FILES:
        src = (S / f).read_text(encoding='utf8')
        toks = lex(src); out = []
        for idx, (k, t) in enumerate(toks):
            if k == 's' and CJK.search(t):
                prev = ''.join(x[1] for x in toks[max(0, idx - 2):idx]).rstrip()
                val = unescape(t[1:-1])
                line = src[:sum(len(x[1]) for x in toks[:idx])].count('\n') + 1
                lits.append({'file': f, 'line': line, 'zh': val})
                if not prev.endswith('tr('):
                    t = 'tr(' + t + ')'
            out.append(t)
        if not dry: (S / f).write_text(''.join(out), encoding='utf8')
    uniq = {}
    for l in lits: uniq.setdefault(l['zh'], []).append(f"{l['file']}:{l['line']}")
    (S / 'i18n' / 'js_literals.json').write_text(json.dumps([{'zh': k, 'where': v[:3]} for k, v in uniq.items()], ensure_ascii=False, indent=1), encoding='utf8')
    print(len(lits), '处字面量,', len(uniq), '条不重复')
main('--dry' in sys.argv)

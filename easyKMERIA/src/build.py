#!/usr/bin/env python3
"""把 src/ 下的各部分合并成单个 static/index.html。python build.py"""
from pathlib import Path
S = Path(__file__).resolve().parent
r = lambda n: (S / n).read_text(encoding='utf8')
import json
(S / 'envsh.js').write_text('/* 由 build.py 从 env-activate.sh 生成, 不要手改 */\nconst ENV_ACTIVATE_SH = ' + json.dumps(r('env-activate.sh').rstrip('\n'), ensure_ascii=False) + ';\n', encoding='utf8')
head = r('head.html').rstrip() + '\n<meta name="easykmera-token" content="__EASYKMERA_TOKEN__">\n'
import base64
def datauri(n): return 'data:image/png;base64,' + base64.b64encode((S / 'assets' / n).read_bytes()).decode()
IMG = {'@@IMG_FLOWER@@': datauri('flower.png'), '@@IMG_ABSTRACT@@': datauri('abstract.png')}
def sub(t):
    for k, v in IMG.items(): t = t.replace(k, v)
    return t
css = sub(r('style.css'))
body = r('nav.html') + '<div class="main-col">\n<main>\n' + r('home.html') + r('gen.body.html') + r('clusters.html') + r('pg.body.html') + r('assistant.html') + r('ai.html') + '</main>\n' + r('foot.html') + '</div>\n'
html = head + '<style>' + css + '</style>\n<div class="app">\n' + sub(body) + '</div>\n'
for tid, f in (('tpl-deploy', 'tpl-deploy.txt'), ('tpl-submit', 'tpl-submit.txt'), ('tpl-postgwas', 'tpl-postgwas.txt'), ('tpl-plotr', 'plot.R.txt')):
    t = r(f)
    assert '</script>' not in t
    html += f'<script type="text/plain" id="{tid}">\n{t.rstrip()}\n</script>\n'
en_path = S / 'i18n' / 'en.json'
i18n_js = r('i18n.js').replace('/*I18N_EN*/{}/*END*/', json.dumps(json.loads(en_path.read_text(encoding='utf8')) if en_path.exists() else {}, ensure_ascii=False))
js = i18n_js + '\n' + r('envsh.js') + '\n' + r('core.js') + '\n' + r('postgwas.js') + '\n(function(){\n' + r('ui-common.js') + r('ui-gen.js') + r('ui-clusters.js') + r('ui-agent.js') + r('ui-pg.js') + '\n})();\n'
html += '<script>\n' + js + '</script>\n'
# 独立运行版: 带完整的 doctype / head
parts = head.split('<style>')[0]
standalone = ('<!doctype html>\n<html lang="zh-CN"><head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n'
              + head + '<style>' + css + '</style>\n</head><body>\n' + html.split('</style>\n', 1)[1] + '</body></html>\n')
(S.parent / 'static' / 'index.html').write_text(standalone, encoding='utf8')
# Artifact 版: 外壳由发布时补上, 这里不写 doctype / head / body
(S.parent / 'dist').mkdir(exist_ok=True)
(S.parent / 'dist' / 'easykmera-preview.html').write_text(html.replace('<meta name="easykmera-token" content="__EASYKMERA_TOKEN__">\n', ''), encoding='utf8')
print('standalone', len(standalone.encode('utf8')) // 1024, 'KB; artifact', len(html.encode('utf8')) // 1024, 'KB')

"""端到端: 真实的 Flask 后端 + 本地 SSH 假服务器(真执行命令) + 假 DeepSeek, 用无头浏览器操作页面。python tests/e2e_ui.py"""
import os, sys, json, socket, subprocess, tempfile, threading, http.server, shutil, time
from pathlib import Path
import paramiko
TMP = Path(tempfile.mkdtemp(prefix='easykmera_e2e_'))
os.environ['EASYKMERA_HOME'] = str(TMP / 'home')
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import app as A

SHOTS = Path(os.environ.get('SHOTS', TMP / 'shots')); SHOTS.mkdir(parents=True, exist_ok=True)
PG = ROOT / 'tests' / 'pg'
# --- 假的集群软件
BIN = TMP / 'bin'; BIN.mkdir()
(BIN / 'sbatch').write_text('#!/bin/bash\necho "sbatch $*" >> ' + str(TMP) + '/sbatch.log\nif [[ "$*" == *--parsable* ]]; then echo $((RANDOM%9000+1000)); else echo "Submitted batch job 4242"; fi\n')
(BIN / 'squeue').write_text("#!/bin/bash\necho '4242|demo|RUNNING|0:10|normal|node01'\n")
(BIN / 'sinfo').write_text("#!/bin/bash\nprintf 'normal*\\nAMD_9654\\n'\n")
for f in BIN.iterdir(): f.chmod(0o755)
ENV = dict(os.environ, PATH=f"{BIN}:{os.environ['PATH']}")
KM = TMP / 'KMERIA'; (KM / 'scripts').mkdir(parents=True); (KM / 'bin').mkdir()
(KM / 'scripts' / 'kmeria_wrapper.pl').write_text(
    '#!/usr/bin/perl\nmy @x=@ARGV; my %a; for(my $i=0;$i<@x;$i++){ if($x[$i] eq "--output"){$a{o}=$x[$i+1]} }\n'
    'foreach my $d (qw(01_kmer_counts 02_kmer_matrices 03_filtered_matrices 04_bimbam 05_association)){ system("mkdir -p $a{o}/$d"); open F,">$a{o}/$d/job_$d.sh"; print F "#!/bin/bash\\necho $d\\n"; close F }\nprint "wrapper ok\\n";\n')
(KM / 'bin' / 'kmeria').write_text('#!/bin/bash\n'); (KM / 'bin' / 'kmeria').chmod(0o755)
FQ = TMP / 'fq'; FQ.mkdir()
for i in range(1, 13):
    for r in ('R1', 'R2'): (FQ / f'SS{i:03d}_{r}.fq.gz').touch()
WORK = TMP / 'work'; WORK.mkdir()

HOST_KEY = paramiko.RSAKey.generate(2048)
class SSH(paramiko.ServerInterface):
    def check_auth_password(self, u, p): return paramiko.AUTH_SUCCESSFUL if (u, p) == ('zhangxingtan', 'secret-pw') else paramiko.AUTH_FAILED
    def get_allowed_auths(self, u): return 'password'
    def check_channel_request(self, k, c): return paramiko.OPEN_SUCCEEDED if k == 'session' else paramiko.OPEN_FAILED_ADMINISTRATIVELY_PROHIBITED
    def check_channel_exec_request(self, ch, cmd):
        threading.Thread(target=run_cmd, args=(ch, cmd.decode()), daemon=True).start(); return True
def run_cmd(ch, cmd):
    data = b''
    while True:
        d = ch.recv(65536)
        if not d: break
        data += d
    p = subprocess.run(['bash', '-c', cmd], input=data, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=ENV)
    ch.sendall(p.stdout); ch.send_exit_status(p.returncode); ch.shutdown(2); ch.close()
def serve(sock):
    while True:
        try: c, _ = sock.accept()
        except OSError: return
        t = paramiko.Transport(c); t.add_server_key(HOST_KEY)
        try: t.start_server(server=SSH())
        except Exception: pass
ss = socket.socket(); ss.bind(('127.0.0.1', 0)); ss.listen(5); SSH_PORT = ss.getsockname()[1]
threading.Thread(target=serve, args=(ss,), daemon=True).start()

class DS(http.server.BaseHTTPRequestHandler):
    def do_POST(self):
        n = int(self.headers.get('Content-Length', 0)); self.rfile.read(n)
        pal = {"options": {"ptSize": "0.4", "alpha": "0.6", "ymax": "9", "labelTop": "3", "stripChr": False, "dpi": "600", "format": "png", "lineP": "5e-7"},
               "palette": {"name": "晨雾蓝绿", "chromosome_colors": ["#264653", "#5FA8A0"], "highlight": "#F4A261", "threshold_line": "#6C757D", "background": "#FBFBF8", "text": "#22333B"},
               "notes": "低饱和蓝绿，高亮暖橙。", "unsupported": ["标注基因名"]}
        self.send_response(200); self.send_header('Content-Type', 'application/json'); self.end_headers()
        self.wfile.write(json.dumps({"choices": [{"message": {"content": json.dumps(pal)}}]}).encode())
    def log_message(self, *a): pass
ds = http.server.HTTPServer(('127.0.0.1', 0), DS); DS_PORT = ds.server_address[1]
threading.Thread(target=ds.serve_forever, daemon=True).start()

PORT = 18765
threading.Thread(target=lambda: A.app.run(host='127.0.0.1', port=PORT, threaded=True), daemon=True).start()
time.sleep(1.2)

from playwright.sync_api import sync_playwright
ok = bad = 0
def check(name, cond, extra=''):
    global ok, bad
    if cond: ok += 1; print('  通过', name)
    else: bad += 1; print('  失败', name, extra)
errs = []
with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={'width': 1360, 'height': 900}, locale='zh-CN')
    pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.on('console', lambda m: errs.append(m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'Failed to load resource' not in m.text else None)
    pg.goto(f'http://127.0.0.1:{PORT}/'); pg.wait_for_timeout(800)
    print('首页'); check('已连接程序', pg.inner_text('#modeChip') == '已连接程序'); check('没有只读提示', not pg.is_visible('#staticNotice'))

    print('集群连接')
    pg.evaluate("location.hash='#clusters'"); pg.wait_for_timeout(300)
    check('空状态提示', '还没有保存的连接' in pg.inner_text('#clusterList'))
    pg.click('#btnAddCluster')
    pg.fill('#cf_name', 'AGIS'); pg.fill('#cf_host', '127.0.0.1'); pg.fill('#cf_port', str(SSH_PORT)); pg.fill('#cf_user', 'zhangxingtan')
    pg.fill('#cf_password', 'Zq9pwTest'); pg.fill('#cf_workdir', str(WORK)); pg.fill('#cf_queue', 'AMD_9654'); pg.fill('#cf_home', str(KM))
    pg.click('#cf_save'); pg.wait_for_selector('.cl-card'); 
    check('卡片显示', 'AGIS' in pg.inner_text('#clusterList') and 'SLURM' in pg.inner_text('#clusterList') and '密码 · 已保存' in pg.inner_text('#clusterList'))
    check('页面里没有明文密码', 'Zq9pwTest' not in pg.content())
    pg.click('[data-a=test]'); pg.wait_for_function("document.querySelector('.cl-res').innerText.length>10 && !document.querySelector('.cl-res').innerText.includes('正在连接')")
    check('错误密码有中文提示', '登录失败' in pg.inner_text('.cl-res'), pg.inner_text('.cl-res'))
    pg.click('[data-a=edit]'); pg.fill('#cf_password', 'secret-pw'); pg.click('#cf_save'); pg.wait_for_timeout(500)
    pg.click('[data-a=test]'); pg.wait_for_function("document.querySelector('.cl-res').innerText.includes('登录成功')")
    t = pg.inner_text('.cl-res')
    check('测试连接成功并识别环境', 'SLURM' in t and 'AMD_9654' in t and '找到了 KMERIA' in t and '工作路径存在' in t, t)
    pg.screenshot(path=str(SHOTS / 'e_clusters.png'), full_page=True)

    print('流程生成并部署到集群')
    pg.evaluate("location.hash='#generate'"); pg.wait_for_timeout(300)
    pg.click('#tabs [data-tab=cluster]'); pg.wait_for_selector('#rg_sel')
    pg.click('[data-act=apply]'); pg.wait_for_timeout(300)
    check('填充了工作路径和安装目录', pg.input_value('#workDir') == str(WORK) and pg.input_value('#kmeriaHome') == str(KM))
    pg.click('#tabs [data-tab=deploy]'); pg.fill('#fastqDir', str(FQ)); pg.wait_for_timeout(400)
    pg.click('#tabs [data-tab=cluster]'); pg.wait_for_selector('#rg_sel')
    pg.click('[data-step="0"]')
    try: pg.wait_for_function("document.getElementById('rg_log').textContent.includes('退出码')", timeout=15000)
    except Exception:
        print('DEBUG html=', pg.evaluate("document.getElementById('p-cluster').innerHTML.slice(0,300)"), 'btn=', pg.evaluate("document.querySelector('[data-step=\"0\"]').disabled"), 'log=', repr(pg.evaluate("document.getElementById('rg_log') && document.getElementById('rg_log').textContent")), 'errs=', errs, 'toast=', pg.inner_text('#toast'))
        raise
    log = pg.inner_text('#rg_log')
    check('部署脚本运行成功', '退出码 0（成功）' in log and '原始数据文件齐全' in log and '预览将要提交的作业' in log, log[-500:])
    proj = WORK / 'ss_TN_demo'
    check('远端生成了输入文件和提交脚本', (proj / 'sample_depth.tsv').exists() and (proj / 'submit_chain.sh').exists())
    pg.click('[data-step="1"]'); check('提交需要二次确认', '再点一次确认' in pg.inner_text('[data-step="1"]'))
    pg.click('[data-step="1"]'); pg.wait_for_function("document.getElementById('rg_log').textContent.includes('已提交')", timeout=60000)
    sb = (TMP / 'sbatch.log').read_text()
    check('按依赖链提交 5 个作业', sb.count('sbatch') == 5 and sb.count('--dependency=afterok:') == 4, sb)
    pg.click('#p-cluster [data-step="2"]'); pg.wait_for_function("document.getElementById('rg_log').textContent.includes('RUNNING')")
    check('查看作业状态', True)
    pg.fill('#rg_dir', str(proj / 'results')); pg.click('[data-act=ls]'); pg.wait_for_selector('#rg_files button[data-ls$="01_kmer_counts"]')
    check('浏览服务器目录', True)
    pg.screenshot(path=str(SHOTS / 'e_generate_cluster.png'))

    print('AI 设置')
    pg.evaluate("location.hash='#ai'"); pg.wait_for_timeout(300)
    pg.click('#btnAiEdit'); pg.fill('#ai_key', 'sk-e2e'); pg.fill('#ai_base', f'http://127.0.0.1:{DS_PORT}')
    pg.screenshot(path=str(SHOTS / 'e_ai_modal.png'))
    pg.click('#ai_save'); pg.wait_for_timeout(500)
    check('已配置', pg.inner_text('#aiChip') == '已配置' and 'sk-e2e' not in pg.content())
    pg.click('#btnAiTest'); pg.wait_for_function("document.getElementById('aiTestMsg').innerText.includes('正常')")
    check('测试连接', True)

    print('Post-GWAS')
    pg.evaluate("location.hash='#postgwas'"); pg.wait_for_timeout(400)
    pg.fill('#pg_assocGlob', str(PG / 'assoc' / '*.assoc.txt')); pg.fill('#pg_genomeFasta', str(PG / 'genome.fa')); pg.fill('#pg_figTitle', 'TN GWAS'); pg.fill('#pg_lineP', '1e-6')
    pg.check('#pg_genesOn'); pg.fill('#pg_gff3', str(PG / 'genes.gff3')); pg.fill('#pg_window', '20000'); pg.wait_for_timeout(400)
    check('检查清单没有错误', '错误' not in pg.inner_text('#pgChecks').split('\n')[0], pg.inner_text('#pgChecks'))
    pg.fill('#pg_aiPrompt', '低饱和蓝绿色，高亮用暖橙，点小一点半透明'); pg.click('#pg_aiGo')
    try: pg.wait_for_function("document.getElementById('pg_aiMsg').innerText.includes('已修改')", timeout=8000)
    except Exception: print('AI 提示内容:', pg.inner_text('#pg_aiMsg'), errs); raise
    check('AI 配色已应用到图', pg.input_value('#pgColors input[data-k=chr0]') == '#264653' and pg.input_value('#pgColors input[data-k=hi]') == '#f4a261')
    pg.screenshot(path=str(SHOTS / 'e_postgwas_ai.png'))
    pg.click('#pgTabs [data-tab=cluster]'); pg.wait_for_selector('#rp_sel')
    pg.click('#pp-cluster [data-step="2"]'); pg.click('#pp-cluster [data-step="2"]')
    pg.wait_for_function("document.getElementById('rp_log').textContent.includes('退出码')", timeout=120000)
    log = pg.inner_text('#rp_log')
    out = proj / 'postgwas'
    check('Post-GWAS 在远端跑通', '退出码 0（成功）' in log and '已画出 21 个' in log, log[-600:])
    check('生成了曼哈顿图(AI 选了只要 PNG, 所以没有 PDF)和基因列表', (out / 'manhattan.png').exists() and not (out / 'manhattan.pdf').exists() and (out / 'sigkmer.near_genes.txt').stat().st_size > 0)
    msg = pg.inner_text('#pg_aiMsg')
    check('AI 改了参数并说明做不到的需求', pg.input_value('#pg_ptSize') == '0.4' and pg.input_value('#pg_alpha') == '0.6' and pg.input_value('#pg_ymax') == '9' and pg.input_value('#pg_labelTop') == '3' and pg.input_value('#pg_lineP') == '5e-07' and not pg.is_checked('#pg_stripChr') and '标注基因名' in msg and '做不到' in msg, msg)
    check('图用了 AI 配色', '"#264653"' in (out / 'plot_manhattan_easy.R').read_text() and '#F4A261' in (out / 'plot_manhattan_easy.R').read_text())
    shutil.copy(out / 'manhattan.png', SHOTS / 'e_remote_manhattan.png')
    pg.click('#pgTabs [data-tab=preview]'); pg.wait_for_timeout(300)
    pg.screenshot(path=str(SHOTS / 'e_postgwas_preview.png'))
    b.close()
print('页面脚本错误:', errs or '无')
print(f'\n通过 {ok} 项, 失败 {bad} 项')
sys.exit(1 if bad or errs else 0)

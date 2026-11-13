"""AI 助手端到端: 真实后端 + 本地 SSH 假服务器 + 假 DeepSeek + 有状态的假 squeue。python tests/e2e_agent.py"""
import os, sys, json, socket, subprocess, tempfile, threading, http.server, time, urllib.request
from pathlib import Path
import paramiko
TMP = Path(tempfile.mkdtemp(prefix='easykmera_ag_'))
os.environ['EASYKMERA_HOME'] = str(TMP / 'home')
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import app as A
SHOTS = Path(os.environ.get('SHOTS', TMP / 'shots')); SHOTS.mkdir(parents=True, exist_ok=True)
PG = ROOT / 'tests' / 'pg'
WORK = TMP / 'work'; WORK.mkdir()
PROJ = WORK / 'agent_demo'
BIN = TMP / 'bin'; BIN.mkdir()
(BIN / 'sbatch').write_text(f'''#!/bin/bash
echo "sbatch $*" >> {TMP}/sbatch.log
if [[ "$*" == *--parsable* ]]; then id=$((RANDOM%9000+1000)); echo $id; else id=4242; echo "Submitted batch job $id"; fi
echo $id >> {TMP}/ids
if [[ "$*" == *postgwas.sh* ]]; then bash postgwas.sh > postgwas.sh.out 2>&1; fi
''')
(BIN / 'squeue').write_text(f'''#!/bin/bash
n=$(grep -c sbatch {TMP}/sbatch.log)
c=$(cat {TMP}/sq_$n 2>/dev/null || echo 0); c=$((c+1)); echo $c > {TMP}/sq_$n
echo 'JOBID|NAME|STATE|TIME|PARTITION|NODELIST(REASON)'
ids=$(tail -5 {TMP}/ids)
if [ $c -le 2 ]; then for i in $ids; do echo "$i|j|RUNNING|0:10|normal|node01"; done
elif [ $n -eq 5 ]; then for i in $ids; do echo "$i|j|PENDING|0:00|normal|(DependencyNeverSatisfied)"; done
elif [ $n -eq 10 ]; then mkdir -p {PROJ}/results/05_association; echo "kmer 1 2" > {PROJ}/results/05_association/out.assoc.txt; fi
''')
(BIN / 'sinfo').write_text("#!/bin/bash\nprintf 'normal*\\nAMD_9654\\n'\n")
(TMP / 'conda_base' / 'etc' / 'profile.d').mkdir(parents=True)
(TMP / 'conda_base' / 'etc' / 'profile.d' / 'conda.sh').write_text('conda() { return 0; }\n')
(BIN / 'conda').write_text(f'#!/bin/bash\nif [ "$1" = info ]; then echo {TMP}/conda_base; fi\nexit 0\n')
(BIN / 'git').write_text('#!/bin/bash\nd="${@: -1}"; mkdir -p "$d/bin" "$d/external_tools" "$d/bimbamAsso" "$d/lib"; printf "#!/bin/bash\\necho kmeria usage: count kctm filter m2b asso\\n" > "$d/bin/kmeria"; touch "$d/kmeria_env.yaml"; echo cloned >> ' + str(TMP) + '/git.log\n')
for f in BIN.iterdir(): f.chmod(0o755)
ENV = dict(os.environ, PATH=f"{BIN}:{os.environ['PATH']}")
KM = TMP / 'KMERIA'; (KM / 'scripts').mkdir(parents=True); (KM / 'bin').mkdir()
(KM / 'scripts' / 'kmeria_wrapper.pl').write_text(
    '#!/usr/bin/perl\nmy @x=@ARGV; my %a; for(my $i=0;$i<@x;$i++){ if($x[$i] eq "--output"){$a{o}=$x[$i+1]} }\n'
    'foreach my $d (qw(01_kmer_counts 02_kmer_matrices 03_filtered_matrices 04_bimbam 05_association)){ system("mkdir -p $a{o}/$d"); open F,">$a{o}/$d/job_$d.sh"; print F "#!/bin/bash\\necho $d\\n"; close F }\nprint "wrapper ok\\n";\n')
(KM / 'bin' / 'kmeria').write_text('#!/bin/bash\n'); (KM / 'bin' / 'kmeria').chmod(0o755)
KM_NEW = TMP / 'apps' / 'KMERIA_new'
FQ = TMP / 'fq'; FQ.mkdir()
for i in range(1, 13):
    for r in ('R1', 'R2'): (FQ / f'SS{i:03d}_{r}.fq.gz').touch()
(TMP / 'ids').write_text(''); (TMP / 'sbatch.log').write_text('')

HOST_KEY = paramiko.RSAKey.generate(2048)
class SSH(paramiko.ServerInterface):
    def check_auth_password(self, u, p): return paramiko.AUTH_SUCCESSFUL if (u, p) == ('lab', 'pw') else paramiko.AUTH_FAILED
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

SEEN = []
class DS(http.server.BaseHTTPRequestHandler):
    def do_POST(self):
        n = int(self.headers.get('Content-Length', 0)); body = json.loads(self.rfile.read(n)); sysmsg = body['messages'][0]['content']; SEEN.append(body)
        if '表单填写助手' in sysmsg:
            out = {"patch": {"project": "agent_demo", "workDir": str(WORK), "kmeriaHome": str(KM), "fastqDir": str(FQ), "ploidy": "4", "scheduler": "slurm", "queue": "AMD_9654", "memory": "64"}, "notes": "按你的描述填了 8 项", "missing": []}
        elif '故障诊断' in sysmsg:
            out = {"summary": "作业因内存不足被杀", "cause": "OOM", "severity": "error", "retry_recommended": True, "fixes": [{"text": "把每个作业的内存调到 256G", "patch": {"memory": "256"}}]}
        else:
            out = {}
        self.send_response(200); self.send_header('Content-Type', 'application/json'); self.end_headers()
        self.wfile.write(json.dumps({"choices": [{"message": {"content": json.dumps(out)}}]}).encode())
    def log_message(self, *a): pass
ds = http.server.HTTPServer(('127.0.0.1', 0), DS); DS_PORT = ds.server_address[1]
threading.Thread(target=ds.serve_forever, daemon=True).start()

PORT = 18766
threading.Thread(target=lambda: A.app.run(host='127.0.0.1', port=PORT, threaded=True), daemon=True).start()
time.sleep(1.2)
def post(path, data):
    r = urllib.request.Request(f'http://127.0.0.1:{PORT}{path}', json.dumps(data).encode(), {'Content-Type': 'application/json', 'X-Token': A.TOKEN}, method='POST')
    return json.loads(urllib.request.urlopen(r).read())
post('/api/clusters', {'name': 'LAB', 'host': '127.0.0.1', 'port': SSH_PORT, 'username': 'lab', 'auth': 'password', 'password': 'pw', 'scheduler': 'slurm', 'queue': 'AMD_9654', 'workdir': str(WORK), 'kmeria_home': str(KM)})
post('/api/ai/settings', {'api_key': 'sk-e2e', 'base_url': f'http://127.0.0.1:{DS_PORT}', 'model': 'deepseek-chat'})

from playwright.sync_api import sync_playwright
ok = bad = 0
def check(name, cond, extra=''):
    global ok, bad
    if cond: ok += 1; print('  通过', name)
    else: bad += 1; print('  失败', name, extra)
errs = []
def ask(pg, text, yes=True, timeout=60000):
    """等到出现含 text 的批准卡片, 点批准或拒绝"""
    sel = f'.ag-item.ask:has-text("{text}")'
    pg.wait_for_selector(sel, timeout=timeout)
    return sel
def click(pg, sel, yes=True): pg.click(f'{sel} [data-a="{1 if yes else 0}"]')

with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1360, 'height': 900}, locale='zh-CN')
    pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.on('console', lambda m: errs.append(m.text) if m.type == 'error' and 'Failed to load resource' not in m.text else None)
    pg.goto(f'http://127.0.0.1:{PORT}/'); pg.evaluate("localStorage.setItem('easykmera.pollMs','400')"); pg.reload(); pg.wait_for_timeout(800)
    pg.evaluate("location.hash='#assistant'"); pg.wait_for_timeout(400)
    check('助手页面可见且选到了连接', pg.is_visible('#v-assistant') and 'LAB' in pg.inner_text('#ag_sel'))

    print('一句话填表')
    pg.fill('#ag_text', '项目 agent_demo，4 倍体，数据在 FQ，slurm 的 AMD_9654 队列'); pg.click('#ag_fill')
    sel = ask(pg, '应用到“流程生成”表单')
    check('展示了现在→改成的对照', 'agent_demo' in pg.inner_text(sel) and '队列' in pg.inner_text(sel))
    check('批准前表单没变', pg.input_value('#project') == 'ss_TN_demo')
    click(pg, sel, False); pg.wait_for_timeout(300)
    check('拒绝后表单仍没变', pg.input_value('#project') == 'ss_TN_demo')
    pg.click('#ag_fill'); sel = ask(pg, '应用到“流程生成”表单'); click(pg, sel); pg.wait_for_selector('.ag-item:has-text("已填入表单")')
    check('批准后表单已填入', pg.input_value('#project') == 'agent_demo' and pg.input_value('#workDir') == str(WORK) and pg.input_value('#queue') == 'AMD_9654' and pg.input_value('#ploidy') == '4')
    sent = json.dumps(SEEN[-1], ensure_ascii=False)
    check('发给 DeepSeek 的内容不含样本表', 'SS001' not in sent and 'S. spontaneum' not in sent)
    pg.screenshot(path=str(SHOTS / 'a_fill.png'), full_page=True)

    print('检查并安装 KMERIA')
    pg.evaluate("location.hash='#generate'"); pg.fill('#kmeriaHome', str(KM_NEW)); pg.wait_for_timeout(300)
    pg.evaluate("location.hash='#assistant'"); pg.click('#ag_clear'); pg.click('#ag_install')
    sel = ask(pg, '安装 KMERIA'); check('安装前没有克隆过', not (TMP / 'git.log').exists())
    check('卡片里展示了完整脚本', 'git clone https://github.com/Sh1ne111/KMERIA.git' in pg.inner_text(sel))
    click(pg, sel); pg.wait_for_selector('.ag-item:has-text("KMERIA 安装完成")', timeout=60000)
    check('安装目录真的建好了', (KM_NEW / 'bin' / 'kmeria').exists() and 'INSTALL_OK' in (WORK / 'easykmera_setup' / 'install_kmeria.sh.log').read_text())
    pg.click('#ag_install'); pg.wait_for_selector('.ag-item:has-text("无需再装")', timeout=30000)
    check('再次检查时识别为已安装', True)
    pg.evaluate("location.hash='#generate'"); pg.fill('#kmeriaHome', str(KM)); pg.wait_for_timeout(300)

    print('一键跑完整流程')
    pg.evaluate("location.hash='#postgwas'"); pg.wait_for_timeout(400)
    pg.fill('#pg_assocGlob', str(PG / 'assoc' / '*.assoc.txt')); pg.fill('#pg_genomeFasta', str(PG / 'genome.fa')); pg.fill('#pg_lineP', '1e-6'); pg.fill('#pg_figTitle', 'AI GWAS'); pg.wait_for_timeout(300)
    pg.evaluate("location.hash='#assistant'"); pg.click('#ag_clear'); pg.click('#ag_auto')
    sel = ask(pg, '部署项目')
    check('批准前没有创建任何东西', not PROJ.exists() and (TMP / 'sbatch.log').read_text() == '')
    click(pg, sel)
    sel = ask(pg, '提交作业')
    check('批准提交前没有提交', (TMP / 'sbatch.log').read_text() == '' and (PROJ / 'submit_chain.sh').exists())
    click(pg, sel)
    sel = ask(pg, '应用修正并重试', timeout=120000)
    diag = pg.inner_text('.ag-item:has-text("诊断结果")')
    check('第一轮失败后给出诊断', '内存不足' in diag and 'OOM' in diag, diag)
    check('提交了 5 个作业', (TMP / 'sbatch.log').read_text().count('sbatch') == 5)
    check('批准前内存仍是 64', pg.input_value('#memory') == '64')
    pg.screenshot(path=str(SHOTS / 'a_retry.png'), full_page=True)
    click(pg, sel); pg.wait_for_timeout(500)
    check('批准后内存改为 256', pg.input_value('#memory') == '256')
    sel = ask(pg, '修正后重新部署'); click(pg, sel)
    sel = ask(pg, '提交作业'); check('重提前提示会重新提交整个流程', '重新提交整个流程' in pg.inner_text(sel)); click(pg, sel)
    sel = ask(pg, '提交 Post-GWAS', timeout=120000)
    check('主流程监控到结果文件', 'out.assoc.txt' in pg.inner_text('#ag_log'))
    check('第二轮共 10 次 sbatch', (TMP / 'sbatch.log').read_text().count('sbatch') == 10)
    click(pg, sel)
    pg.wait_for_selector('.ag-item:has-text("结果已就绪")', timeout=240000)
    out = PROJ / 'postgwas'
    check('Post-GWAS 在远端产出曼哈顿图', (out / 'manhattan.png').exists() and (out / 'sigkmer.assoc.blastn.plot.txt').stat().st_size > 0)
    check('结果已填入 Post-GWAS 预览', pg.evaluate("document.getElementById('pg_plotText').value.length") > 20 and pg.evaluate("document.querySelector('input[name=pg_src]:checked').value") == 'mine')
    pg.screenshot(path=str(SHOTS / 'a_done.png'), full_page=True)
    pg.evaluate("location.hash='#postgwas'"); pg.wait_for_timeout(500); pg.screenshot(path=str(SHOTS / 'a_pg.png'))

    print('停止按钮')
    pg.evaluate("location.hash='#assistant'"); pg.click('#ag_clear'); pg.click('#ag_auto')
    sel = ask(pg, '部署项目'); check('任务中出现停止按钮且其他按钮被禁用', pg.is_visible('#ag_stop') and pg.is_disabled('#ag_auto'))
    pg.click('#ag_stop'); pg.wait_for_timeout(500)
    check('停止后卡片标为已拒绝并恢复按钮', '已拒绝' in pg.inner_text(sel.replace('.ask', '')) and not pg.is_visible('#ag_stop') and not pg.is_disabled('#ag_auto'))
    b.close()
print('页面脚本错误:', errs or '无')
print(f'\n通过 {ok} 项, 失败 {bad} 项')
sys.exit(1 if bad or errs else 0)

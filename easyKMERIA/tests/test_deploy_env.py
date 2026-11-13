"""部署脚本的环境激活: 模拟 micromamba 只在被中断的 .bashrc 里出现、作业在干净环境里运行。python tests/test_deploy_env.py"""
import json, os, subprocess, sys, tempfile
from pathlib import Path
S = Path(__file__).resolve().parent.parent / 'src'
T = Path(tempfile.mkdtemp(prefix='easykmeria_env_'))
HOME = T / 'home'; HOME.mkdir()
ok = bad = 0
def check(n, c, x=''):
    global ok, bad
    if c: ok += 1; print('  通过', n)
    else: bad += 1; print('  失败', n, x)

# 假环境: micromamba 在非常规目录, 只有 .bashrc 里写了路径, 且 .bashrc 在交互判断处就 return 了
MM = HOME / 'opt/mm/bin/micromamba'; MM.parent.mkdir(parents=True)
MM.write_text('#!/bin/bash\nif [ "$1 $2" = "shell hook" ]; then cat <<\'X\'\nmicromamba() { if [ "$1" = activate ]; then export CONDA_DEFAULT_ENV=$2; export PATH="$MAMBA_ROOT_PREFIX/envs/$2/bin:$PATH"; fi; }\nX\nfi\n'); MM.chmod(0o755)
KB = HOME / 'opt/mm/root/envs/kmeriaenv/bin'; KB.mkdir(parents=True)
(KB / 'kmc').write_text('#!/bin/bash\necho "kmc ok"\n'); (KB / 'kmc').chmod(0o755)
(HOME / '.bashrc').write_text(f"bind -x '\"\\C-l\":x' 2>/dev/null\n[[ $- != *i* ]] && return\nexport MAMBA_EXE='{MM}'\nexport MAMBA_ROOT_PREFIX='{HOME}/opt/mm/root'\n")
KM = T / 'KMERIA'; (KM / 'scripts').mkdir(parents=True); (KM / 'bin').mkdir()
(KM / 'bin/kmeria').write_text('#!/bin/bash\n'); (KM / 'bin/kmeria').chmod(0o755)
(KM / 'scripts/kmeria_wrapper.pl').write_text(
    '#!/usr/bin/perl\nmy @x=@ARGV; my %a; for(my $i=0;$i<@x;$i++){ if($x[$i] eq "--output"){$a{o}=$x[$i+1]} }\n'
    'foreach my $d (qw(01_kmer_counts 02_kmer_matrices)){ system("mkdir -p $a{o}/$d"); open F,">$a{o}/$d/job_$d.sh"; print F "#!/bin/bash\\n#SBATCH -J x\\n#SBATCH -p normal\\n\\nkmc --version\\nkmeria\\n"; close F }\nprint "wrapper ok\\n";\n')
FQ = T / 'fq'; FQ.mkdir()
for i in range(1, 13):
    for r in ('R1', 'R2'): (FQ / f'SS{i:03d}_{r}.fq.gz').touch()
js = 'const tr=x=>x, trTpl=x=>x;\n' + (S / 'envsh.js').read_text() + (S / 'core.js').read_text() + (S / 'postgwas.js').read_text() + f"""
const cfg = Object.assign({{}}, EXAMPLE, {{workDir:{json.dumps(str(T/'work'))}, kmeriaHome:{json.dumps(str(KM))}, fastqDir:{json.dumps(str(FQ))}, condaEnv:'kmeriaenv'}});
const r = generate(cfg, {{deploy: {json.dumps((S/'tpl-deploy.txt').read_text())}, submit: {json.dumps((S/'tpl-submit.txt').read_text())}}});
process.stdout.write(r.deploy);"""
deploy = subprocess.run(['node', '-e', js], capture_output=True, text=True).stdout
(T / 'deploy.sh').write_text(deploy)
check('生成了脚本且语法正确', subprocess.run(['bash', '-n', str(T / 'deploy.sh')]).returncode == 0)
env = dict(PATH='/usr/bin:/bin', HOME=str(HOME), LANG='C.UTF-8')
r = subprocess.run(['bash', str(T / 'deploy.sh')], capture_output=True, text=True, env=env)
out = r.stdout + r.stderr
check('部署脚本成功', r.returncode == 0, out[-600:])
check('从 .bashrc 里写的路径找到 micromamba 并激活', '[环境] 已用 micromamba 激活 kmeriaenv' in out, out[-600:])
proj = next((T / 'work').iterdir())
check('作业脚本被加上 source 且 #SBATCH 头仍在最前', all(
    (lambda L: L[1].startswith('#SBATCH') and 'easykmeria_env.sh' in '\n'.join(L) and L.index(next(x for x in L if 'easykmeria_env.sh' in x)) > 2)((proj / 'results' / d / f'job_{d}.sh').read_text().splitlines())
    for d in ('01_kmer_counts', '02_kmer_matrices')))
check('提示处理了 2 个作业脚本', '已给 2 个作业脚本加上环境加载' in out, out[-400:])
# 模拟计算节点: 干净环境, 只有系统 PATH
job = subprocess.run(['bash', str(proj / 'results/01_kmer_counts/job_01_kmer_counts.sh')], capture_output=True, text=True, env=dict(PATH='/usr/bin:/bin', HOME=str(HOME)))
check('作业在干净环境里也找得到 kmc 和 kmeria', 'kmc ok' in job.stdout and 'command not found' not in job.stderr, job.stdout + job.stderr)
r2 = subprocess.run(['bash', str(T / 'deploy.sh')], capture_output=True, text=True, env=env)
txt = (proj / 'results/01_kmer_counts/job_01_kmer_counts.sh').read_text()
check('重复部署不会重复添加', r2.returncode == 0 and txt.count('easykmeria_env.sh') == 1, txt)
# 作业脚本带 set -eu 时也不能被环境文件弄坏
(proj / 'results/01_kmer_counts/strict.sh').write_text('#!/bin/bash\nset -euo pipefail\nsource "%s/easykmeria_env.sh"\nkmc --version\n' % proj)
j2 = subprocess.run(['bash', str(proj / 'results/01_kmer_counts/strict.sh')], capture_output=True, text=True, env=dict(PATH='/usr/bin:/bin', HOME=str(HOME)))
check('严格模式 (set -eu) 的作业也能加载', j2.returncode == 0 and 'kmc ok' in j2.stdout, j2.stdout + j2.stderr)
# 什么都找不到时: 只警告, 部署照常完成
HOME2 = T / 'home2'; HOME2.mkdir()
r3 = subprocess.run(['bash', str(T / 'deploy.sh')], capture_output=True, text=True, env=dict(PATH='/usr/bin:/bin', HOME=str(HOME2)))
check('找不到环境管理器时给出指引但不中断', '没有找到 micromamba / mamba / conda' in r3.stdout + r3.stderr and 'type -a micromamba conda' in r3.stdout + r3.stderr and r3.returncode == 0, (r3.stdout + r3.stderr)[-500:])
print(f'\n通过 {ok} 项, 失败 {bad} 项'); sys.exit(1 if bad else 0)

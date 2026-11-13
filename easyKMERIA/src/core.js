
const DEFAULTS = {
  project:'', workDir:'', kmeriaHome:'', condaEnv:'kmeriaenv', preCmd:'',
  dataMode:'dir', fastqDir:'', sampleText:'', defaultDepth:'30', readLen:'150', depthIsSeq:true,
  phenoText:'', trait:'0',
  ploidy:'4', kmer:'31', counter:'kmeria', assoTool:'bimbam', minAbund:'4', maxAbund:'1000', missing:'0.6', step:'all',
  batchSize:'4', kmcMemory:'16', kctmBatch:'10000', sketchSize:'8000000', bgzfThreads:'16', kinPrec:'10', outPrec:'5',
  covar:'', kinship:'', sepStrand:true,
  scheduler:'slurm', queue:'share', threads:'32', memory:'64', time:'720:00:00'
};
const EXAMPLE = Object.assign({}, DEFAULTS, {
  project:'ss_TN_demo', workDir:'/share/home/yourname/kmeria_runs', kmeriaHome:'/share/apps/KMERIA',
  fastqDir:'/share/data/rawdata/S_spontaneum',
  sampleText:['sample\tdepth','SS001\t31.2','SS002\t28.7','SS003\t35.4','SS004\t26.9','SS005\t30.1','SS006\t33.8',
    'SS007\t29.5','SS008\t27.3','SS009\t34.6','SS010\t31.9','SS011\t25.8','SS012\t32.4'].join('\n'),
  phenoText:['sample\tTN\tAP','SS001\t18\t83.2','SS002\t25\t79.5','SS003\t12\t88.1','SS004\t31\t74.6','SS005\t22\t81.0',
    'SS006\t15\t85.7','SS007\t28\t77.3','SS008\t20\t80.4','SS009\t9\t90.2','SS010\t26\t76.8','SS011\tNA\t82.5','SS012\t17\tNA'].join('\n'),
  ploidy:'8', queue:'normal'
});

const NA_TOKENS = new Set(['','na','nan','n/a','-','.','null']);
const STEP_DIR = {count:'01_kmer_counts', kctm:'02_kmer_matrices', filter:'03_filtered_matrices', m2b:'04_bimbam', asso:'05_association'};
const STEP_ORDER = ['count','kctm','filter','m2b','asso'];

function isNum(s){ return s !== '' && s != null && !isNaN(Number(s)); }
function cleanLines(text){
  return String(text || '').replace(/\r/g,'').split('\n').map(l => l.replace(/\s+$/,'')).filter(l => l.trim() !== '' && !l.trim().startsWith('#'));
}
function detectDelim(lines){ const j = lines.join('\n'); if (j.includes('\t')) return '\t'; if (j.includes(',')) return ','; return ' '; }
function splitLine(line, d){
  if (d === '\t') return line.split('\t').map(s => s.trim());
  if (d === ',') return line.split(',').map(s => s.trim());
  return line.trim().split(/\s+/);
}
function joinPath(a, b){ return String(a || '').replace(/\/+$/,'') + '/' + String(b || '').replace(/^\/+/,''); }
function fmtNum(x){ return String(Number(Number(x).toFixed(2))); }

function parseSampleTable(text, mode){
  const lines = cleanLines(text);
  const d = detectDelim(lines);
  const entries = []; let header = false;
  lines.forEach((ln, i) => {
    const c = splitLine(ln, d);
    if (i === 0) {
      const first = c[0] || '';
      if (/^(sample|samples|sample_?id|sample_?name|id|name|样本|样本名|样本名称)$/i.test(first) ||
          (mode === 'dir' && c.length > 1 && c[1] !== '' && !isNum(c[1]))) { header = true; return; }
    }
    if (mode === 'perSample') entries.push({name:c[0]||'', r1:c[1]||'', r2:c[2]||'', depth:c[3]||'', ploidy:c[4]||''});
    else entries.push({name:c[0]||'', r1:'', r2:'', depth:c[1]||'', ploidy:c[2]||''});
  });
  return {entries: entries.filter(e => e.name !== ''), header};
}

function parsePheno(text){
  const lines = cleanLines(text);
  if (!lines.length) return {traits:[], rows:new Map(), header:false};
  const d = detectDelim(lines);
  const cells = lines.map(l => splitLine(l, d));
  const first = cells[0];
  const header = first.length > 1 && first.slice(1).some(v => !NA_TOKENS.has(v.toLowerCase()) && !isNum(v));
  const data = header ? cells.slice(1) : cells;
  const ncol = Math.max(0, ...cells.map(c => c.length)) - 1;
  const traits = [];
  for (let i = 0; i < ncol; i++) traits.push(header && first[i+1] ? first[i+1] : 'trait' + (i+1));
  const rows = new Map();
  data.forEach(c => {
    if (!c[0]) return;
    const vals = [];
    for (let i = 0; i < ncol; i++) {
      const v = (c[i+1] == null ? '' : c[i+1]);
      vals.push(NA_TOKENS.has(v.toLowerCase()) ? 'NA' : v);
    }
    rows.set(c[0], vals);
  });
  return {traits, rows, header};
}

function fill(tpl, map){ return tpl.replace(/@@([A-Z_]+)@@/g, (m, k) => (k in map ? map[k] : m)); }

function buildWrapper(cfg, sh, ctx){
  const V = sh ? {HOME:'$KMERIA_HOME', PROJ:'$PROJECT_DIR', OUT:'$OUT_DIR', FASTQ:'$FASTQ_DIR'}
               : {HOME:cfg.kmeriaHome.replace(/\/+$/,''), PROJ:ctx.projectDir, OUT:ctx.outDir, FASTQ:ctx.fastqDir};
  const q = s => sh ? '"' + s + '"' : s;
  const step = cfg.step;
  const useKmc = cfg.counter === 'kmc';
  const a = [];
  const add = (flag, val) => a.push(val === undefined ? flag : flag + ' ' + val);
  const opt = (flag, val, def) => { if (String(val) !== String(def) && String(val) !== '') add(flag, val); };

  // 官方 wrapper 把 --output 当作"基础目录", 自己在下面建 01_kmer_counts ... 05_association (单步模式也一样),
  // 所以 --output 永远只给 results, 不要给 results/05_association (否则会变成 05_association/05_association)。
  const inOut = {
    all:   [q(V.FASTQ), q(V.OUT)],
    count: [q(V.FASTQ), q(V.OUT)],
    kctm:  [q(V.OUT + '/' + STEP_DIR.count), q(V.OUT)],
    filter:[q(V.OUT + '/' + STEP_DIR.kctm), q(V.OUT)],
    m2b:   [q(V.OUT + '/' + STEP_DIR.filter), q(V.OUT)],
    asso:  [q(V.OUT + '/' + STEP_DIR.m2b), q(V.OUT)]
  }[step];
  add('--step', step);
  add('--input', inOut[0]);
  add('--output', inOut[1]);
  add('--samples', q(V.PROJ + '/sample.list'));
  add('--threads', cfg.threads);
  add('--memory', cfg.memory + 'G');
  add('--time', cfg.time);
  add('--scheduler', cfg.scheduler);
  if (cfg.scheduler !== 'local' && cfg.queue) add('--queue', cfg.queue);

  const doCount = step === 'all' || step === 'count';
  const doKctm = step === 'all' || step === 'kctm';
  const doFilter = step === 'all' || step === 'filter';
  const doM2b = step === 'all' || step === 'm2b';
  const doAsso = step === 'all' || step === 'asso';
  if (doCount) {
    add('--kmer', cfg.kmer);
    add('--min-abund', cfg.minAbund);
    add('--batch-size', cfg.batchSize);
    if (useKmc) { add('--use-kmc'); add('--kmc-memory', cfg.kmcMemory); }
    else if (cfg.sepStrand) add('-C');
  }
  if (doCount || doFilter) add('--max-abund', cfg.maxAbund);
  if (doKctm) opt('--kctm-batch', cfg.kctmBatch, DEFAULTS.kctmBatch);
  if (doFilter) { add('--ploidy', cfg.ploidy); add('--missing', cfg.missing); }
  if (doFilter || doM2b) add('--depth-file', q(V.PROJ + '/sample_depth.tsv'));
  if (doM2b) { opt('--bgzf-threads', cfg.bgzfThreads, DEFAULTS.bgzfThreads); opt('--sketch-size', cfg.sketchSize, DEFAULTS.sketchSize); }
  if (doAsso) {
    add('--pheno', q(V.PROJ + '/phenotypes.txt'));
    add('--pheno-col', String(Number(cfg.trait || 0) + 1));
    if (cfg.covar) add('--covar', q(cfg.covar));
    if (cfg.kinship) add('--kinship', q(cfg.kinship));
    if (cfg.assoTool === 'bimbam') {
      add('--use-bimbam-tools');
      opt('--kinship-precision', cfg.kinPrec, DEFAULTS.kinPrec);
      opt('--output-precision', cfg.outPrec, DEFAULTS.outPrec);
    }
  }
  const head = 'perl ' + q(V.HOME + '/scripts/kmeria_wrapper.pl');
  return head + ' \\\n' + a.map(x => '    ' + x).join(' \\\n');
}

function generate(cfg, T){
  cfg = Object.assign({}, DEFAULTS, cfg);
  const checks = [];
  const add = (level, msg) => checks.push({level, msg});
  const step = cfg.step;
  const needFastq = step === 'all' || step === 'count';
  const needPheno = step === 'all' || step === 'asso';

  const projectDir = joinPath(cfg.workDir || '/path/to/workdir', cfg.project || 'kmeria_project');
  const outDir = projectDir + '/results';
  const fastqDir = cfg.dataMode === 'dir' ? (cfg.fastqDir || '/path/to/fastq').replace(/\/+$/,'') : projectDir + '/00_fastq_links';
  const ctx = {projectDir, outDir, fastqDir};

  // ---- 样本 ----
  const parsed = parseSampleTable(cfg.sampleText, cfg.dataMode);
  const seen = new Set(), entries = [], dups = [];
  parsed.entries.forEach(e => { if (seen.has(e.name)) dups.push(e.name); else { seen.add(e.name); entries.push(e); } });
  const pheno = parsePheno(cfg.phenoText);
  let names = entries.map(e => e.name);
  if (!names.length && pheno.rows.size) { names = Array.from(pheno.rows.keys()); }

  if (!cfg.project) add('err', tr('请填写项目名。'));
  else if (/[^A-Za-z0-9_.-]/.test(cfg.project)) add('warn', tr('项目名里有特殊字符，建议只用字母、数字、下划线。'));
  if (!cfg.workDir) add('err', tr('请填写工作目录。'));
  if (!cfg.kmeriaHome) add('err', tr('请填写 KMERIA 安装路径。'));
  const pathFields = [[tr('工作目录'), cfg.workDir], [tr('KMERIA 安装路径'), cfg.kmeriaHome], [tr('原始数据目录'), cfg.dataMode === 'dir' ? cfg.fastqDir : ''], [tr('协变量文件'), cfg.covar], [tr('亲缘矩阵文件'), cfg.kinship]];
  pathFields.forEach(([n, p]) => {
    if (!p) return;
    if (/\s/.test(p)) add('err', n + tr('里有空格，请换成没有空格的路径。'));
    else if (!/^[\/~$]/.test(p)) add('warn', n + tr('不是绝对路径，建议以 / 开头。'));
  });
  if (!names.length) add('err', tr('还没有识别到样本，请在第 2 步填写样本表。'));
  else add('ok', tr('识别到 ') + names.length + tr(' 个样本。'));
  if (dups.length) add('err', tr('样本名重复：') + dups.slice(0,4).join('、') + (dups.length > 4 ? tr(' 等') : '') + tr('，已只保留第一次出现的。'));
  const bad = names.filter(n => /[^A-Za-z0-9._-]/.test(n));
  if (bad.length) add('warn', tr('样本名含特殊字符（') + bad.slice(0,3).join('、') + tr('），建议只用字母、数字、下划线、点和短横线。'));
  if (names.length && names.length < 30) add('warn', tr('样本数只有 ') + names.length + tr(' 个。GWAS 通常需要上百个样本，样本太少统计效能很低。'));

  if (needFastq && cfg.dataMode === 'dir' && !cfg.fastqDir) add('err', tr('请填写原始数据目录。'));
  if (cfg.dataMode === 'perSample') {
    const miss = entries.filter(e => !e.r1);
    if (miss.length) add('err', tr('有 ') + miss.length + tr(' 个样本没有 R1 路径（如 ') + miss[0].name + tr('）。每行格式：样本名 R1路径 R2路径。'));
    const noR2 = entries.filter(e => e.r1 && !e.r2);
    if (noR2.length) add('warn', noR2.length + tr(' 个样本没有 R2 路径，将按单个文件处理。官方文档只说明了双端数据的命名方式。'));
  }

  // ---- 深度与倍性 ----
  const L = Number(cfg.readLen) || 150, k = Number(cfg.kmer) || 31;
  const depthLines = [];
  const hasPloidyCol = entries.some(e => isNum(e.ploidy));
  let noDepth = 0, badDepth = 0;
  names.forEach(n => {
    const e = entries.find(x => x.name === n) || {};
    let d = isNum(e.depth) ? Number(e.depth) : NaN;
    if (isNaN(d)) { d = Number(cfg.defaultDepth); noDepth++; }
    if (!isFinite(d) || d <= 0) { badDepth++; d = 0; }
    const kd = cfg.depthIsSeq ? d * (L - k + 1) / L : d;
    const cols = [n, fmtNum(kd)];
    if (hasPloidyCol) cols.push(isNum(e.ploidy) ? String(Number(e.ploidy)) : String(cfg.ploidy));
    depthLines.push(cols.join('\t'));
  });
  if (badDepth) add('err', badDepth + tr(' 个样本没有可用的测序深度，请在样本表或「默认测序深度」里填写。'));
  else if (names.length && noDepth === names.length) add('warn', tr('所有样本都使用默认测序深度 ') + cfg.defaultDepth + tr('X，请确认。'));
  else if (noDepth) add('warn', noDepth + tr(' 个样本没有写深度，使用默认值 ') + cfg.defaultDepth + 'X。');
  if (cfg.depthIsSeq && names.length) add('ok', tr('k-mer 深度 = 测序深度 × ') + fmtNum((L - k + 1) / L) + tr('（reads ') + L + tr(' bp，k=') + k + tr('）。'));
  if (hasPloidyCol) add('ok', tr('样本表带有倍性列，深度文件按混合倍性格式写出（三列）。'));
  if (k > 31) add('warn', tr('k-mer 长度 ') + k + tr(' 超过 31，官方不建议用于 GWAS。'));
  if (cfg.counter === 'kmc') add('ok', tr('使用 KMC 计数，需要 kmc 和 kmc_tools 在 PATH 里（通常在 external_tools 目录）。'));
  if (cfg.scheduler === 'local') add('warn', tr('选择了直接运行。不要在登录节点上跑，请先申请计算节点。'));
  else if (!cfg.queue) add('warn', tr('没有填写队列，将使用 wrapper 的默认队列 share。'));

  // ---- 表型 ----
  const traitIdx = Math.max(0, Math.min(Number(cfg.trait || 0), Math.max(0, pheno.traits.length - 1)));
  cfg.trait = String(traitIdx);
  const phenoLines = [];
  let nPheno = 0;
  names.forEach(n => {
    const v = pheno.rows.get(n);
    if (v) nPheno++;
    const vals = pheno.traits.map((_, i) => v ? v[i] : 'NA');
    phenoLines.push([n].concat(vals.length ? vals : ['NA']).join('\t'));
  });
  if (!pheno.rows.size) {
    if (needPheno) add('err', tr('请在第 3 步填写表型表。'));
  } else {
    const missingP = names.length - nPheno;
    if (missingP > 0) add('warn', missingP + tr(' 个样本在表型表里找不到，已填 NA。'));
    const extra = Array.from(pheno.rows.keys()).filter(n => !seen.has(n) && names.length && !names.includes(n));
    if (extra.length) add('warn', tr('表型表里有 ') + extra.length + tr(' 个样本不在样本表中，已忽略（如 ') + extra[0] + tr('）。'));
    const nNA = names.filter(n => { const v = pheno.rows.get(n); return v && v[traitIdx] === 'NA'; }).length;
    if (nNA) add('warn', tr('性状 ') + pheno.traits[traitIdx] + tr(' 有 ') + nNA + tr(' 个样本是 NA。'));
    if (nPheno) add('ok', tr('表型已对应 ') + nPheno + tr(' 个样本，当前分析性状：') + pheno.traits[traitIdx] + tr('（第 ') + (traitIdx + 1) + tr(' 列）。'));
  }

  // ---- 脚本 ----
  const date = new Date().toISOString().slice(0, 10);
  const sampleList = names.join('\n');
  const depthFile = depthLines.join('\n');
  const phenoFile = phenoLines.join('\n');

  const envBlock = buildEnvBlock(cfg.condaEnv, cfg.preCmd);

  let linkBlock = '';
  if (cfg.dataMode === 'perSample') {
    const extOf = p => { const m = p.match(/\.(f(?:ast)?[aq])(\.gz)?$/i); return m ? '.' + m[1].toLowerCase() + (m[2] || '') : '.fq.gz'; };
    const L2 = ['', tr('# ---------- 2b. 为每个样本建立规范命名的数据链接 ----------'),
      tr('# KMERIA 只认 样本名_R1 / 样本名_R2 这样的文件名, 这里用软链接统一命名, 不会复制或改动原始数据'),
      'mkdir -p "$FASTQ_DIR"'];
    entries.forEach(e => {
      if (!e.r1) return;
      L2.push('ln -sfn "' + e.r1 + '" "$FASTQ_DIR/' + e.name + '_R1' + extOf(e.r1) + '"');
      if (e.r2) L2.push('ln -sfn "' + e.r2 + '" "$FASTQ_DIR/' + e.name + '_R2' + extOf(e.r2) + '"');
    });
    L2.push(tr('echo ">>> 已建立数据链接: $FASTQ_DIR"'), '');
    linkBlock = L2.join('\n');
  }

  const extraCmds = [];
  if (cfg.counter === 'kmc' && (step === 'all' || step === 'count')) extraCmds.push('kmc', 'kmc_tools');
  if (step === 'all' || step === 'filter' || step === 'm2b' || step === 'asso') extraCmds.push('plink');
  if (cfg.assoTool === 'gemma' && (step === 'all' || step === 'asso')) extraCmds.push('gemma');

  let dataCheck = '';
  if (needFastq) {
    dataCheck = [
      tr('echo ">>> 检查原始数据是否齐全"'),
      'missing_n=0',
      'while IFS= read -r s; do',
      '  [ -z "$s" ] && continue',
      '  if ! compgen -G "$FASTQ_DIR/${s}_R1.f*q*" >/dev/null && ! compgen -G "$FASTQ_DIR/${s}_1.f*q*" >/dev/null; then',
      tr('    echo "  [缺失] 样本 $s 在 $FASTQ_DIR 里没有找到 ${s}_R1.* 或 ${s}_1.* 文件"'),
      '    missing_n=$((missing_n + 1))',
      '  fi',
      'done < sample.list',
      tr('if [ "$missing_n" -gt 0 ]; then fail=1; else echo "  原始数据文件齐全"; fi')
    ].join('\n');
  }

  const stepDirs = (step === 'all' ? STEP_ORDER : [step]).map(s => STEP_DIR[s]);
  const submit = fill(T.submit.replace(/^\n/, '').replace(/\n$/, ''), {OUT_DIR: outDir, SCHED: cfg.scheduler, STEP_DIRS: stepDirs.join(' ')});
  const wrapperSh = buildWrapper(cfg, true, ctx);
  const wrapperLit = buildWrapper(cfg, false, ctx);
  const deploy = fill(T.deploy.replace(/^\n/, '').replace(/\n$/, ''), {
    PROJECT: cfg.project || tr('(未命名)'), DATE: date, KMERIA_HOME: cfg.kmeriaHome.replace(/\/+$/,'') || '/path/to/KMERIA',
    PROJECT_DIR: projectDir, FASTQ_DIR: cfg.dataMode === 'dir' ? fastqDir : '$PROJECT_DIR/00_fastq_links',
    ENV_BLOCK: envBlock, JOB_ENV: buildJobEnv(cfg).replace(/\n$/, ''), SAMPLE_LIST: sampleList, DEPTH_FILE: depthFile, PHENO_FILE: phenoFile,
    LINK_BLOCK: linkBlock, CHECK_DATA_LABEL: needFastq ? tr('和原始数据') : '', EXTRA_CMDS: extraCmds.join(' ') || 'true',
    DATA_CHECK: dataCheck, WRAPPER_CMD: wrapperSh, SUBMIT_SCRIPT: submit
  });

  const home = cfg.kmeriaHome.replace(/\/+$/,'') || '/path/to/KMERIA';
  const install = [
    tr('# 1. 下载 KMERIA 到你填写的安装路径'),
    'git clone https://github.com/Sh1ne111/KMERIA.git "' + home + '"',
    'cd "' + home + '"',
    '',
    tr('# 2. 创建 conda 环境 (有 mamba 的话把 conda 换成 mamba 会快很多)'),
    'conda env create -f kmeria_env.yaml',
    'conda activate kmeriaenv',
    '',
    tr('# 3. 设置库路径, 赋予可执行权限, 加入 PATH'),
    'export LD_LIBRARY_PATH="' + home + '/lib:${LD_LIBRARY_PATH:-}"',
    'chmod 755 "' + home + '"/bin/* "' + home + '"/external_tools/* "' + home + '"/bimbamAsso/*',
    'export PATH="' + home + '/bin:' + home + '/bimbamAsso:' + home + '/external_tools:$PATH"',
    '',
    tr('# 4. 验证: 能看到 count / kctm / filter / m2b / asso 等命令列表就说明装好了'),
    'kmeria'
  ].join('\n');

  const traitNames = pheno.traits;
  return {checks, files:{sampleList, depthFile, phenoFile}, wrapperLit, deploy, submit, install,
          summary:{n:names.length, nPheno, traits:traitNames, traitIdx, header:parsed.header}};
}

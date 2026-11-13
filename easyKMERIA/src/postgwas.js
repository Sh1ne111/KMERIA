/* ===== Post-GWAS: 显著 k-mer 定位与曼哈顿图 ===== */
const PALETTES = [
  {name:tr('海洋蓝'), chr:['#2F5D8A','#7FB3D5'], hi:'#E4572E', line:'#6B7280', bg:'#FFFFFF', text:'#1F2933'},
  {name:tr('Nature 蓝紫'), chr:['#3C5488','#8491B4'], hi:'#E64B35', line:'#6B7280', bg:'#FFFFFF', text:'#1F2933'},
  {name:tr('森林绿'), chr:['#2E6F4E','#8CBF9F'], hi:'#D1495B', line:'#6B7280', bg:'#FFFFFF', text:'#1F2933'},
  {name:tr('石板灰'), chr:['#4A5568','#A0AEC0'], hi:'#C05621', line:'#718096', bg:'#FFFFFF', text:'#1A202C'},
  {name:tr('暖赭色'), chr:['#8C5E3C','#D9A66B'], hi:'#1B6CA8', line:'#6B7280', bg:'#FFFFFF', text:'#2D2A26'},
  {name:tr('深色背景'), chr:['#5FB3B3','#9AA9B8'], hi:'#FFB347', line:'#8A98A6', bg:'#12181B', text:'#E6EDF0'}
];

const PG_DEFAULTS = {
  assocGlob:'', thresh:'1e-4', lineP:'1e-4',
  dbMode:'genome', genomeFasta:'', blastDb:'', evalue:'1e-3', wordSize:'9', blastThreads:'8',
  excludeRe:'Chr00|tig|scaf', useFai:true,
  genesOn:false, gff3:'', window:'100000',
  figW:'12', figH:'4.5', figTitle:'',
  ptSize:'0.7', hiSize:'1', alpha:'1', ymax:'', showLine:true, lineStyle:'dashed', line2P:'', grid:false, stripChr:true,
  shape:'circle', dpi:'300', format:'both', labelTop:'0', fontScale:'1',
  pgMemory:'16', pgTime:'24:00:00',
  palette: PALETTES[0]
};

const SAFE_PATH = /^[A-Za-z0-9_.\/~+@%=:,*?\[\]{}-]+$/;
const SAFE_RE = /^[A-Za-z0-9_|.^$()*+?\[\]-]+$/;
const SAFE_NUM = /^(?:\d+(?:\.\d+)?|\.\d+)(?:[eE][-+]?\d+)?$/;
const HEX6 = /^#[0-9a-fA-F]{6}$/;

const ENV_CALL = env => '_easy_u=0; case "$-" in *u*) _easy_u=1; set +u ;; esac\n_easy_activate "' + env + '" || true\nif [ "$_easy_u" = 1 ]; then set -u; fi\n';
function buildEnvBlock(condaEnv, preCmd){
  let s = '';
  if ((preCmd || '').trim()) s += tr('# 额外环境命令\n') + preCmd.trim() + '\n\n';
  const env = (condaEnv || '').trim();
  if (env) s += trTpl(ENV_ACTIVATE_SH) + '\n' + ENV_CALL(env);
  return s;
}
/* 作业脚本在计算节点上 source 的环境文件: 额外命令 + 激活环境 + KMERIA 的库路径和 PATH */
function buildJobEnv(cfg){
  const home = (cfg.kmeriaHome || '').replace(/\/+$/, '') || '/path/to/KMERIA';
  let s = tr('# easyKMERIA 生成: 每个作业脚本开头会 source 这个文件, 让计算节点上也能找到 KMERIA 和它的依赖 (kmc 等)\n');
  s += buildEnvBlock(cfg.condaEnv, cfg.preCmd);
  s += 'export LD_LIBRARY_PATH="' + home + '/lib:${LD_LIBRARY_PATH:-}"\n';
  s += 'export PATH="' + home + '/bin:' + home + '/bimbamAsso:' + home + '/external_tools:$PATH"\n';
  return s;
}

function schedHeader(sched, queue, threads, memG, time, name){
  const L = [];
  if (sched === 'slurm') {
    L.push('#SBATCH -J ' + name);
    if (queue) L.push('#SBATCH -p ' + queue);
    L.push('#SBATCH -N 1', '#SBATCH -n 1', '#SBATCH --cpus-per-task=' + threads, '#SBATCH --mem=' + memG + 'G',
      '#SBATCH --time=' + time, '#SBATCH -o ' + name + '.%j.log', '#SBATCH -e ' + name + '.%j.err');
  } else if (sched === 'sge') {
    L.push('#$ -N ' + name, '#$ -cwd');
    if (queue) L.push('#$ -q ' + queue);
    L.push('#$ -pe smp ' + threads, '#$ -o ' + name + '.log', '#$ -e ' + name + '.err');
  } else if (sched === 'pbs') {
    L.push('#PBS -N ' + name);
    if (queue) L.push('#PBS -q ' + queue);
    L.push('#PBS -l nodes=1:ppn=' + threads, '#PBS -l walltime=' + time, '#PBS -o ' + name + '.log', '#PBS -e ' + name + '.err');
  }
  return L.join('\n');
}

function rQuote(s){ return String(s || '').replace(/["\\\r\n]/g, ''); }

function generatePostGwas(pg, gen, T){
  pg = Object.assign({}, PG_DEFAULTS, pg);
  gen = Object.assign({scheduler:'slurm', queue:'', condaEnv:'', preCmd:'', assoTool:'bimbam', project:'', workDir:''}, gen);
  const checks = [];
  const add = (level, msg) => checks.push({level, msg});
  const projectDir = joinPath(gen.workDir || '/path/to/workdir', gen.project || 'kmeria_project');
  const outDir = projectDir + '/postgwas';

  if (!gen.project || !gen.workDir) add('err', tr('请先在「流程生成」里填写项目名和工作目录，Post-GWAS 的文件会放在项目目录的 postgwas 子目录里。'));
  const glob = (pg.assocGlob || '').trim();
  if (!glob) add('err', tr('请填写关联分析结果文件（可以用通配符，例如 ') + projectDir + tr('/results/05_association/*.assoc.txt）。'));
  else if (!SAFE_PATH.test(glob)) add('err', tr('关联结果路径里有空格或特殊字符，请换成只含字母、数字、下划线、点、短横线和 / 的路径。'));
  else add('warn', tr('请先到 05_association 目录用 ls 核对结果文件的实际名字，通配符匹配不到文件时脚本会直接报错退出。'));

  const thresh = (pg.thresh || '').trim();
  if (!SAFE_NUM.test(thresh) || Number(thresh) <= 0 || Number(thresh) >= 1) add('err', tr('p 值阈值要写成 0 到 1 之间的数，例如 1e-4。'));
  const lineP = (pg.lineP || '').trim() || thresh;
  if (!SAFE_NUM.test(lineP) || Number(lineP) <= 0 || Number(lineP) >= 1) add('err', tr('曼哈顿图的阈值线要写成 0 到 1 之间的数。'));
  else if (Number(lineP) > Number(thresh)) add('warn', tr('阈值线比筛选阈值更宽松，图上所有点都会被当作显著。'));

  const genome = (pg.genomeFasta || '').trim(), db = (pg.blastDb || '').trim();
  if (pg.dbMode === 'genome') {
    if (!genome) add('err', tr('请填写参考基因组 fasta 的路径，脚本会自动建 BLAST 库。'));
    else if (!SAFE_PATH.test(genome)) add('err', tr('参考基因组路径里有空格或特殊字符。'));
  } else {
    if (!db) add('err', tr('请填写已有 BLAST 库的路径（不带 .nhr 等后缀）。'));
    else if (!SAFE_PATH.test(db)) add('err', tr('BLAST 库路径里有空格或特殊字符。'));
  }
  const excl = (pg.excludeRe || '').trim();
  if (excl && !SAFE_RE.test(excl)) add('err', tr('排除未定位序列的表达式只能含字母、数字和 | . ^ $ ( ) * + ? [ ] -'));
  if (!SAFE_NUM.test((pg.evalue || '').trim())) add('err', tr('BLAST e-value 要写成数字，例如 1e-3。'));
  if (!/^\d+$/.test(String(pg.wordSize).trim()) || Number(pg.wordSize) < 4) add('err', tr('word_size 要写成不小于 4 的整数。'));
  if (!/^\d+$/.test(String(pg.blastThreads).trim()) || Number(pg.blastThreads) < 1) add('err', tr('BLAST 线程数要写成正整数。'));
  if (pg.genesOn) {
    if (!(pg.gff3 || '').trim()) add('err', tr('要查找附近的基因，请填写 GFF3 注释文件路径。'));
    else if (!SAFE_PATH.test(pg.gff3.trim())) add('err', tr('GFF3 路径里有空格或特殊字符。'));
    if (!/^\d+$/.test(String(pg.window).trim())) add('err', tr('基因搜索范围要写成整数（单位 bp）。'));
  }
  const pal = pg.palette || PALETTES[0];
  if (![pal.hi, pal.line, pal.bg, pal.text].every(c => HEX6.test(c)) || !(pal.chr || []).length || !pal.chr.every(c => HEX6.test(c))) add('err', tr('配色里有不是 #RRGGBB 格式的颜色。'));
  const title = rQuote(pg.figTitle);
  if (/[^\x00-\x7F]/.test(title)) add('warn', tr('图标题含有中文，R 默认字体可能显示成方框，建议用英文标题。'));
  if (!(Number(pg.figW) > 0 && Number(pg.figH) > 0)) add('err', tr('图片宽高要写成正数。'));
  const rng = (v, lo, hi) => SAFE_NUM.test(String(v).trim()) && Number(v) >= lo && Number(v) <= hi;
  if (!rng(pg.ptSize, 0.1, 5) || !rng(pg.hiSize, 0.1, 8)) add('err', tr('点大小要写成数字（普通点 0.1–5，高亮点 0.1–8）。'));
  if (!rng(pg.alpha, 0.05, 1)) add('err', tr('点的不透明度要在 0.05 到 1 之间。'));
  if (String(pg.ymax).trim() && !rng(pg.ymax, 1, 1000)) add('err', tr('纵轴上限要写成数字，留空表示自动。'));
  if (String(pg.line2P).trim() && !(SAFE_NUM.test(String(pg.line2P).trim()) && Number(pg.line2P) > 0 && Number(pg.line2P) < 1)) add('err', tr('第二条阈值线要写成 0 到 1 之间的 p 值，留空表示不画。'));
  if (!/^\d+$/.test(String(pg.dpi).trim()) || Number(pg.dpi) < 72 || Number(pg.dpi) > 1200) add('err', tr('分辨率（dpi）要写成 72 到 1200 之间的整数。'));
  if (!/^\d+$/.test(String(pg.labelTop).trim()) || Number(pg.labelTop) > 30) add('err', tr('标注的 k-mer 个数要写成 0 到 30 的整数。'));
  if (!rng(pg.fontScale, 0.5, 3)) add('err', tr('文字缩放要写成 0.5 到 3 之间的数。'));
  if (gen.scheduler !== 'local' && !gen.queue) add('warn', tr('没有填写队列，提交时需要你自己在脚本头部补上队列。'));
  if (gen.assoTool === 'gemma') add('ok', tr('按 GEMMA 的结果格式取列（k-mer 在第 2 列，p 值在第 12 列）。'));
  else add('ok', tr('按 bimbamAsso 的结果格式取列（k-mer 在第 1 列，p 值在第 3 列）。'));

  const gem = gen.assoTool === 'gemma';
  const kCol = gem ? 2 : 1, pCol = gem ? 12 : 3;
  const useFai = pg.dbMode === 'genome' && pg.useFai && genome;
  const extra = [];
  if (pg.dbMode === 'genome') extra.push('makeblastdb');
  if (useFai) extra.push('samtools');
  if (pg.genesOn) extra.push('bedtools');

  let dbBlock;
  if (pg.dbMode === 'genome') {
    dbBlock = tr('DB="$OUT_DIR/ref_db/ref"\nif ! ls "${DB}".n* >/dev/null 2>&1; then\n  echo ">>> 用 makeblastdb 建库: ') + genome + '"\n  mkdir -p ref_db\n  makeblastdb -in "' + genome + '" -dbtype nucl -out "$DB"\nfi';
  } else {
    dbBlock = 'DB="' + db + tr('"\nls "${DB}".n* >/dev/null 2>&1 || { echo "[缺失] 找不到 BLAST 库: $DB"; exit 1; }');
  }
  const filterBlock = excl
    ? tr('# 去掉未定位的 contig/scaffold, 每个 k-mer 只保留第一个比对结果\n{ grep -vE \'') + excl + '\' sigkmer.assoc.blastn || true; } | awk \'!a[$1]++\' > sigkmer.assoc.blastn.filter'
    : tr('# 每个 k-mer 只保留第一个比对结果\nawk \'!a[$1]++\' sigkmer.assoc.blastn > sigkmer.assoc.blastn.filter');
  const faiBlock = useFai ? tr('\n# 用基因组各条染色体的长度来画满整个横轴\n[ -f "') + genome + '.fai" ] || samtools faidx "' + genome + '"\n' : '';
  const geneBlock = pg.genesOn ? [
    '', tr('# ---------- 5b. 定位 k-mer 附近的基因 (easyKMERIA 扩展步骤, 不属于官方 Wiki 流程) ----------'),
    tr('# 要求 GFF3 里的染色体名和参考基因组一致'),
    'awk \'BEGIN{FS="\\t"; OFS="\\t"} $3=="gene"{id=$9; if (match($9,/ID=[^;]+/)) id=substr($9,RSTART+3,RLENGTH-3); print $1,$4-1,$5,id}\' "' + pg.gff3.trim() + '" | sort -k1,1 -k2,2n > genes.bed',
    'awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $2,$3-1,$3,$1,$4}\' sigkmer.assoc.blastn.plot.txt | sort -k1,1 -k2,2n > sigkmer.bed',
    'bedtools window -a sigkmer.bed -b genes.bed -w ' + pg.window.trim() + ' > sigkmer.near_genes.txt',
    tr('echo ">>> 附近有基因的记录数: $(wc -l < sigkmer.near_genes.txt)"'), ''].join('\n') : '';
  const geneEcho = pg.genesOn ? tr('\necho "  sigkmer.near_genes.txt          显著 k-mer 附近 ±') + pg.window.trim() + tr(' bp 内的基因"') : '';

  const threads = String(pg.blastThreads).trim();
  const script = fill(T.postgwas.replace(/^\n/, '').replace(/\n$/, ''), {
    HEADER: schedHeader(gen.scheduler, gen.queue, threads, pg.pgMemory, pg.pgTime, 'easykmera_postgwas'),
    PROJECT: gen.project || tr('(未命名)'), DATE: new Date().toISOString().slice(0, 10),
    OUT_DIR: outDir, SCRIPT_DIR: '$OUT_DIR', ENV_BLOCK: buildEnvBlock(gen.condaEnv, gen.preCmd) ? '\n' + buildEnvBlock(gen.condaEnv, gen.preCmd) : '',
    EXTRA_CMDS: extra.length ? ' ' + extra.join(' ') : '', ASSOC_GLOB: glob || '/path/to/assoc_results',
    P_COL: String(pCol), K_COL: String(kCol), THRESH: thresh, ASSOC_FILE: 'asso_thres' + thresh + '.txt',
    COL_NOTE: gem ? tr('GEMMA 结果: k-mer 在第 2 列, p 值在第 12 列') : tr('bimbamAsso 结果: k-mer 在第 1 列, p 值在第 3 列'),
    DB_BLOCK: dbBlock, THREADS: threads, EVALUE: (pg.evalue || '').trim(), WORDSIZE: String(pg.wordSize).trim(),
    FILTER_BLOCK: filterBlock, FAI_BLOCK: faiBlock, GENE_BLOCK: geneBlock, GENE_ECHO: geneEcho,
    LINE_P: lineP, FAI_ARG: useFai ? '"' + genome + '.fai"' : 'NONE', EXCL_ARG: excl ? "'" + excl + "'" : 'NONE'
  });
  const rscript = fill(T.plotR.replace(/\n$/, ''), {
    LINE_P: lineP, CHR_COLORS: (pal.chr || []).map(c => '"' + c + '"').join(', '), HI: pal.hi, LINE: pal.line, BG: pal.bg, TEXT: pal.text,
    W: String(Number(pg.figW) || 12), H: String(Number(pg.figH) || 4.5), TITLE: title,
    PT_CEX: String(Number(pg.ptSize) || 0.7), HI_CEX: String(Number(pg.hiSize) || 1), ALPHA: String(Number(pg.alpha) || 1),
    YMAX: String(pg.ymax).trim() && Number(pg.ymax) > 0 ? String(Number(pg.ymax)) : 'NA', SHOW_LINE: pg.showLine ? 'TRUE' : 'FALSE',
    LTY: {solid:'1', dashed:'2', dotted:'3'}[pg.lineStyle] || '2', LINE_TWO: String(pg.line2P).trim() && Number(pg.line2P) > 0 ? String(Number(pg.line2P)) : 'NA',
    GRID: pg.grid ? 'TRUE' : 'FALSE', STRIP: pg.stripChr ? 'TRUE' : 'FALSE', PCH: {circle:'19', triangle:'17', square:'15', hollow:'1'}[pg.shape] || '19',
    DPI: String(parseInt(pg.dpi, 10) || 300), FMT: ['both', 'pdf', 'png'].includes(pg.format) ? pg.format : 'both',
    LABEL_TOP: String(parseInt(pg.labelTop, 10) || 0), FONT: String(Number(pg.fontScale) || 1)
  });
  return {checks, script, rscript, outDir, files:[{name:'postgwas.sh', content:script}, {name:'plot_manhattan_easy.R', content:rscript}]};
}

/* ---- 曼哈顿图预览: 与 R 脚本相同的排版逻辑 ---- */
function natOrder(names){
  const key = n => { const d = String(n).replace(/[^0-9]/g, ''); return d === '' ? NaN : Number(d); };
  return names.slice().sort((a, b) => {
    const ka = key(a), kb = key(b), na = isNaN(ka), nb = isNaN(kb);
    if (na !== nb) return na ? 1 : -1;
    if (!na && ka !== kb) return ka - kb;
    return a < b ? -1 : a > b ? 1 : 0;
  });
}
function parsePlotTable(text){
  const rows = [];
  cleanLines(text).forEach(l => {
    const c = l.split('\t').length >= 4 ? l.split('\t') : l.trim().split(/\s+/);
    if (c.length < 4) return;
    const pos = Number(c[2]), p = Number(c[3]);
    if (!isFinite(pos) || !isFinite(p)) return;
    rows.push({id:c[0], chr:c[1], pos, p: p <= 0 ? 5e-324 : p});
  });
  return rows;
}
function layoutManhattan(rows, faiText, excl, lineP){
  let chrs;
  if (faiText && faiText.trim()) {
    let f = cleanLines(faiText).map(l => l.split(/\s+/)).filter(c => c.length >= 2 && isFinite(Number(c[1]))).map(c => ({chr:c[0], len:Number(c[1])}));
    if (excl) { try { const re = new RegExp(excl); f = f.filter(x => !re.test(x.chr)); } catch (e) {} }
    chrs = f;
  } else {
    const mx = new Map();
    rows.forEach(r => mx.set(r.chr, Math.max(mx.get(r.chr) || 0, r.pos)));
    chrs = Array.from(mx, ([chr, len]) => ({chr, len}));
  }
  const order = natOrder(chrs.map(c => c.chr));
  chrs = order.map(n => chrs.find(c => c.chr === n));
  let off = 0;
  chrs.forEach(c => { c.off = off; off += c.len; });
  const idx = new Map(chrs.map((c, i) => [c.chr, i]));
  const thr = Number(lineP) || 1e-4;
  const pts = [];
  rows.forEach(r => {
    if (!idx.has(r.chr)) return;
    const i = idx.get(r.chr);
    pts.push({x: r.pos + chrs[i].off, y: -Math.log10(r.p), ci: i, sig: r.p < thr, id: r.id, chr: r.chr, pos: r.pos, p: r.p});
  });
  return {chrs, pts, total: off, thr};
}
function demoPlotRows(){
  let s = 20260611;
  const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
  const rows = [], fai = []; let n = 0;
  for (let c = 1; c <= 20; c++) {
    const len = 60e6 + Math.floor(rnd() * 90e6);
    fai.push('Chr' + c + '\t' + len);
    const m = 8 + Math.floor(rnd() * 18);
    for (let i = 0; i < m; i++) {
      const p = Math.pow(10, -(4 + Math.pow(rnd(), 3) * 3.2));
      rows.push({id:'asso_kmer' + (++n), chr:'Chr' + c, pos:Math.floor(rnd() * len), p});
    }
    if (c === 7 || c === 13) {
      const centre = len * (c === 7 ? .42 : .7);
      for (let i = 0; i < 26; i++) rows.push({id:'asso_kmer' + (++n), chr:'Chr' + c, pos:Math.floor(centre + (rnd() - .5) * 1.5e6), p:Math.pow(10, -(5.5 + rnd() * 4.5))});
    }
  }
  return {rows, fai: fai.join('\n')};
}

const fs=require('fs'),path=require('path');
global.tr=x=>x;global.trTpl=x=>x;
const S=path.join(__dirname,'../../src');
const code=fs.readFileSync(S+'/envsh.js','utf8')+fs.readFileSync(S+'/core.js','utf8')+fs.readFileSync(S+'/postgwas.js','utf8');
const m=new Function(code+';return {generatePostGwas,PG_DEFAULTS,PALETTES,layoutManhattan,parsePlotTable,demoPlotRows,natOrder};')();
const T={postgwas:fs.readFileSync(S+'/tpl-postgwas.txt','utf8'),plotR:fs.readFileSync(S+"/plot.R.txt","utf8")};module.exports={m,T};
if(require.main===module){
  const base=process.argv[2], mode=process.argv[3]||'genome';
  const pg=Object.assign({},m.PG_DEFAULTS,{assocGlob:base+'/assoc/*.assoc.txt',genomeFasta:base+'/genome.fa',genesOn:true,gff3:base+'/genes.gff3',window:'20000',lineP:'1e-6',figTitle:'Demo trait'});
  if(mode==='db'){pg.dbMode='db';pg.blastDb=base+'/mydb/ref';pg.useFai=false}
  const r=m.generatePostGwas(pg,{scheduler:'local',project:'p1',workDir:base+'/w',condaEnv:'',assoTool:'bimbam'},T);
  console.log(r.checks.map(c=>c.level+': '+c.msg).join('\n'));
  fs.mkdirSync(r.outDir,{recursive:true});
  r.files.forEach(f=>fs.writeFileSync(r.outDir+'/'+f.name,f.content));
}

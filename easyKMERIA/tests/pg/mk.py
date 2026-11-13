import random, os
random.seed(7)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))
chrs = {'Chr1':400000,'Chr2':300000,'Chr10':250000,'Chr00':50000,'tig001':30000}
g = {k:seq(v) for k,v in chrs.items()}
with open('genome.fa','w') as f:
    for k,s in g.items():
        f.write('>'+k+'\n')
        for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')
# GFF3
with open('genes.gff3','w') as f:
    f.write('##gff-version 3\n')
    f.write('Chr2\tsrc\tgene\t150001\t152000\t.\t+\t.\tID=GeneA;Name=a\n')
    f.write('Chr2\tsrc\tmRNA\t150001\t152000\t.\t+\t.\tID=GeneA.1;Parent=GeneA\n')
    f.write('Chr10\tsrc\tgene\t90001\t93000\t.\t-\t.\tID=GeneB\n')
comp = {'A':'T','C':'G','G':'C','T':'A'}
rc = lambda s: ''.join(comp[c] for c in reversed(s))
rows = []
def add(chrom, pos, p, rev=False):
    k = g[chrom][pos:pos+31]
    if rev: k = rc(k)
    rows.append((k,p))
# 显著: Chr2 峰, Chr10 峰, Chr00/tig 应被排除
for i in range(12): add('Chr2', 150500+i*40, 10**-(6+random.random()*3), rev=(i%3==0))
for i in range(8): add('Chr10', 91000+i*100, 10**-(5+random.random()*2))
add('Chr1', 1000, 5e-5)
add('Chr00', 100, 1e-7); add('tig001', 200, 1e-7)
# 非显著背景
for i in range(300): add('Chr1', random.randint(0,399000), random.uniform(0.001,1))
# 无法比对的随机 k-mer 显著
rows.append((seq(31), 1e-8))
random.shuffle(rows)
os.makedirs('assoc', exist_ok=True)
half = len(rows)//2
for name, part in (('a.assoc.txt', rows[:half]), ('b.assoc.txt', rows[half:])):
    with open('assoc/'+name,'w') as f:
        f.write('kmer\tbeta\tp_wald\n')   # 表头
        for k,p in part: f.write(f'{k}\t0.1\t{p:.3e}\n')
print(len(rows), 'rows')

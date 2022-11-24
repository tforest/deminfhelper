contigs = ["CM026847.1","JADDRP010000578.1","JADDRP010000058.1","CM026856.1"]
line="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AP021	ATA127	JF5447	JF5455	JF5488b	JF5538"

print(contigs)
print(line)

POP="POP:"+(','.join(line.split("\t")[9:]))
contigs=' '.join(contigs)

import subprocess

def smcpp (contigs,POP):
    subprocess.run(['/home/tforest/work/birdsdemography/smcp.sh','contigs','POP'],shell=True)

smcpp (contigs,POP)

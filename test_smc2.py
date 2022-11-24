contigs = ["CM026847.1","JADDRP010000578.1","JADDRP010000058.1","CM026856.1"]
line="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AP021	ATA127	JF5447	JF5455	JF5488b	JF5538"

print(contigs)
print(line)

POP="POP:"+(','.join(line.split("\t")[9:]))
#contigs=' '.join(contigs)

import subprocess
import os
#qsub -cwd -V -N SMCPP -o smcpp.out -e smcpp.err -q short.q

#smc++ vcf2smc /home/tforest/work/birdsdemography/data/Contemporary_Merged.vcf.gz /home/tforest/work/birdsdemography/data/chr1.smc.gz $VARIABLE $2 ;

QSUB_PARAMS = {'JOB_NAME':'SMCPP', 'STD_OUT':'smcpp.out', 'STD_ERR':'smcpp.err', 'QUEUE':'short.q'}
#"qsub -cwd -V -N SMCPP -o smcpp.out -e smcpp.err -q short.q"



def smcpp (contigs,POP):
    submitted_jobs_ids = []
    #print(contigs)
    #subprocess.run(['/home/tforest/work/birdsdemography/smcp.sh','contigs','POP'],shell=True)
    for contig in contigs:
        QSUB_PARAMS['JOB_NAME'] = "SMCPP_"+contig
        QSUB_PARAMS['STD_OUT'] = "smcpp_"+contig+".out"
        QSUB_PARAMS['STD_ERR'] = "smcpp_"+contig+".err"
        #os.system("/usr/local/genome/Anaconda3/condabin/conda init bash")
        cmd = "qsub -cwd -V -N " + QSUB_PARAMS['JOB_NAME'] + " -o " + QSUB_PARAMS['STD_OUT']+ " -e " + QSUB_PARAMS['STD_ERR'] + " -q " + \
        QSUB_PARAMS['QUEUE'] + " -b y '/home/tforest/work/.conda/envs/smcpp/bin/smc++ vcf2smc /home/tforest/work/birdsdemography/data/Contemporary_Merged.vcf.gz /home/tforest/work/birdsdemography/data/"+contig+".smc.gz " + contig + " " + POP+"'"
        output = subprocess.check_output([cmd], shell=True)
        job_submit_output = output.decode(encoding="utf-8")
        submitted_jobs_ids = job_submit_output.split(" ")[2] #liste de l'ID des jobs
	cmd2 = "qsub -cwd -V -N " + QSUB_PARAMS['JOB_NAME']+"_inf" + " -o " + QSUB_PARAMS['STD_OUT']+"_inf" + " -e " + QSUB_PARAMS['STD_ERR']+"_inf" + " -q " + \
        QSUB_PARAMS['QUEUE'] + "-hold-jid " + submitted_jobs_ids " -b y '/home/tforest/work/.conda/envs/smcpp/bin/smc++ smc++ estimate -o "/home/tforest/work/birdsdemography/data/test_smcpp/"+contig+"model.final.json" 4.2e-9 "/home/tforest/work/birdsdemography/data/"+contig+".smc.gz" + \
	print(cmd2)
        #os.system(QSUB_PARAMS+" smc++ vcf2smc /home/tforest/work/birdsdemography/data/Contemporary_Merged.vcf.gz /home/tforest/work/birdsdemography/data/chr1.smc.gz " \
        #+ contig + " " + POP)
    print(submitted_jobs_ids)
smcpp (contigs,POP)

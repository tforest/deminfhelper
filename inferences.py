#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Creates the input files to run the inferences softwares and run the softwares
"""

import re
import os
import subprocess
import itertools
import multiprocessing
import dadi

def input_dadi(popid, sfs, folded, n, out_dir, mask = True):
    #writes the input file to run dadi
    #n est l'effectif de la pop
    #sfs est une np.array avec chaque colonne correspondant à un effectif de variant et chaque ligne
    #correspondant à une population
    #mask = True si masked, False sinon ()
    #folded = True si folded, False sinon
    n=n+1
    #Première ligne :
    dim =str(n)+' '+["unfolded","folded"][folded]
    #ligne de sfs :
    sfs_write=list(itertools.chain.from_iterable([["fs[",str(sfs[i]),"] "] for i in range (len(sfs))]))
    sfs_write=(''.join(map(str,sfs_write))) #concaténation de la liste en une chaîne de caractères
    sfs_write=sfs_write[0:len(sfs_write)-1]#retrait du dernier espace de la ligne
    #écriture du fichier :
    if folded == False:
        with open(out_dir+popid+'_dadi_input.txt', 'w') as f:
            f.write(dim+'\n'+sfs_write+"\n"+["0","1"][mask])
    if folded == True:
        with open(out_dir+popid+'_dadi_input.txt', 'w') as f:
            f.write(dim+'\n'+sfs_write+" fs[0]"*(n-1)+"\n"+["0","1"][mask])


from numpy import array
from matplotlib import pyplot as plt

def dadi_inf(popid,out_dir,out_dir_d,dict,p,lower_bound,upper_bound,p0,mu,L,gen):

    sfs_folded_hir = dict

    print("SFS recent=", sfs_folded_hir)

    sfs_list = list(sfs_folded_hir.values())
    fs = dadi.Spectrum(sfs_list, mask_corners = False)

    fs.to_file(out_dir_d+"SFS_output.fs")

    # folding SFS (applying a mask )
    fs = fs.fold()
    ns = fs.sample_sizes

    pts_l = [ns[0]+5,ns[0]+15,ns[0]+25]

    #Model Definition : Three epoch
    func = dadi.Demographics1D.three_epoch

    param_names= ("nuB","nuF","TB","TF")

    #default : lower_bound = [1, 1, 0.05, 0.01]
    #lower_bound = [0.2, 1, 0.05, 0.01]

    #upper_bound = [10, 4, 0.1, 10]
    #default : p0 = [0.01,0.001, 0.01, 0.01]
    #p0 = [1,1, 0.01, 0.01]

    #   1. Make extrapolation function:
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    #   2. Perturb parameters:
    p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
    #   3. Optimize:

    maxiter = 5000

    popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                    lower_bound=lower_bound,
                                    upper_bound=upper_bound,
                                    verbose=len(p0), maxiter=maxiter, output_file = out_dir_d+"output_"+p+".dadi")

    print("POPT",popt)
    #  4. Calculate best fit model:
    model = func_ex(popt, ns, pts_l)
    #  5. Compute likelihood of the data given the model allele frequency spectrum:
    ll_model = dadi.Inference.ll_multinom(model, fs)
    #  6. Compute the likelihood of the data to itself (best possible LL):
    ll_data=dadi.Inference.ll_multinom(fs, fs)
    #  7. Calculate the best fit theta:
    theta = dadi.Inference.optimal_sfs_scaling(model, fs)
    #  8. Model specific scaling of parameters (will depend on mu and L that you supply):
    Nanc=theta / (4*float(mu)*float(L))
    nu_scaled_dip=popt[0]*Nanc
    print(nu_scaled_dip)
    T_scaled_gen=popt[1]*2*Nanc
    scaled_param_names=("Nanc_FromTheta_scaled_dip","nu_scaled_dip","T_scaled_gen")
    scaled_popt=(Nanc,nu_scaled_dip,T_scaled_gen*float(gen))
    fT = open(out_dir_d+"popt_"+popid+"_dadi.txt", 'w')
    fT.writelines(str(scaled_popt))
    fT.close()

    print("scaled popt", scaled_popt)
    print("theta", theta)
    print("T_scaled_gen", scaled_popt[2])

def run_dadi_cli(popid, out_dir, sfs_path, optimizations=1000, lower_bounds = [0.1, 0.1, 0.1, 0.1],
                 upper_bounds = [50, 5, 30, 10]):
    # create dadi file
    cmd1 = "".join(["dadi-cli InferDM --fs ", sfs_path, \
                    " --model three_epoch --lbounds ", " ".join(map(str, lower_bounds)), \
                    " --ubounds ", " ".join(map(str, upper_bounds)), " --output ", \
                    out_dir+str(popid)+".dadi", " --optimizations ", str(optimizations)])
    print(cmd1)
    os.system(cmd1)
    

def input_stairwayplot2(popid, nseq, L, whether_folded, SFS, mu, year_per_generation, stairway_plot_dir, output_path, temp_blueprint):
    #writes the input file to run stairwayplot2
    #popid: name of the population
    #nseq: number of sequences (2*n for diploids)
    #L: total number of sites?
    #whether_folded: if the vcf is phased or not
    #SFS: sfs as a list without the monomorphic sites
    #mu: mutation rate
    #year_per_generation: generation time
    #stairway_plot_dir: path to stairwayplot2
    locals()['project_dir'] = output_path+popid #where the output of stairwayplot2 will be
    locals()['plot_title'] = popid #name of the plots output by stairwayplot2
    with open(temp_blueprint, "r") as temp, open(output_path+str(popid)+".blueprint","w") as out_file:
        line = temp.readline()
        while line != '':
            if line.split(':')[0] in locals().keys():
                if line.split(':')[0] == 'SFS':
                    out_file.write('sfs: ' + ' '.join(map(str, SFS)) + '\n')
                else:
                    out_file.write(line.split(':')[0]+': '+ str(locals().get(line.split(':')[0]))+'\n')
            elif line.split(':')[0] == "nrand" :
                out_file.write('nrand: '+str(int((nseq-2)/4))+' '+ str(int((nseq-2)/2))+' '+ str(int((nseq-2)*3/4))+' '+ str(int(nseq-2))+'\n')
            else:
                out_file.write(line)
            line = temp.readline()

def run_stairwayplot2(popid, out_dir, path_to_stairwayplot2):
    cmd1 = "".join(["java -cp ", path_to_stairwayplot2, " Stairbuilder ", out_dir, popid, ".blueprint"])
    output = subprocess.check_output([cmd1], shell=True)
    output = output.decode(encoding="utf-8")
    print(output)
    cmd2 = "".join([ "bash ", out_dir, popid, ".blueprint.sh"])
    os.system(cmd2)


def msmc(n, line, header):
    if header :
        if ("##fileformat" in line) :
            for k in range(n):
                vcf_msmc = open("VCF_input_MSMC_"+str(k)+".txt",mode='w')
                vcf_msmc.write(line)
                vcf_msmc.close()
        if ("##FILTER" in line) :
            for k in range(n):
                vcf_msmc = open("VCF_input_MSMC_"+str(k)+".txt",mode='a')
                vcf_msmc.write(line)
                vcf_msmc.close()
        if ("#CHROM" in line):
            header ='\t'.join(line.split("\t")[:9])
            indiv='\t'.join(line[:-1].split("\t")[9:])
            for k in range(n):
                vcf_msmc = open("VCF_input_MSMC_"+str(k)+".txt",mode='a')
                vcf_msmc.write(header+'\t'+indiv.split("\t")[k]+'\n')
                vcf_msmc.close()
    else :
        gens=re.findall('./.',line)
        for k in range(n):
            vcf_msmc = open("VCF_input_MSMC_"+str(k)+".txt",mode='a')
            vcf_msmc.write('\t'.join(line.split("\t")[0:9])+"\t"+gens[k]+"\n")
            vcf_msmc.close()

def msmc2(contigs, popid, pop_ind, vcf, out_dir, mu, gen_time, num_cpus=4):
    #TODO
    #contigs = contigs[:4]
    if len(contigs) == 0:
        print("Error! No contigs to use! Make sure the threshold matches your data.")

    # Get the available CPUs
    available_cpus = list(range(multiprocessing.cpu_count()))
    # Choose the first 'num_cpus' CPUs
    cpu_affinity = available_cpus[:num_cpus]

    for contig in contigs:
        # split in separated vcfs
        cmd2 = " ".join(["bcftools", "view", "-t", contig, vcf, "|",
                "bcftools", "query", "-", "-f", "'%INFO/DP\n'", "|",
                "awk '{ sum += $1 } END { print sum/NR }'"])
        print(cmd2)
        with open(out_dir+contig+"_mean_DP.txt", 'w') as log:
            process = subprocess.Popen(cmd2, shell=True, stdout=log)
            # Set CPU affinity
            process_pid = process.pid
            os.sched_setaffinity(process_pid, cpu_affinity)
    process.wait()
    # need to add asynchronous. Some jobs finish before others
    # IMPORTANT : Causes crashes with empty files and cannot convert string to float
    for contig in contigs:
        # now that we have coverage, run in parallel the vcf.gz splitting
        with open(out_dir+contig+"_mean_DP.txt", 'r') as filin:
            meanDP=float(filin.readline().strip())
        minDP = meanDP/2
        maxDP = meanDP*2
        # omit sites with missing data "./."
        cmd3 = " ".join(["bcftools view -g ^miss -t", contig, vcf,
        "|vcftools --vcf - --minDP", str(minDP), "--maxDP", str(maxDP),
        "--recode --stdout | gzip -c >", out_dir+contig+".vcf.gz"])
        print(cmd3)
        process = subprocess.Popen(cmd3, shell=True)
        # Set CPU affinity
        process_pid = process.pid
        os.sched_setaffinity(process_pid, cpu_affinity)
    process.wait()
    deminfhelper_directory = os.path.dirname(os.path.abspath(__file__))
    for contig in contigs:
        cmd4 = " ".join(["python3",
        deminfhelper_directory+"/scripts/msmc-tools/generate_multihetsep.py",
        contig+".vcf.gz", ">", out_dir+contig+"_msmc_input.txt"])
        print(cmd4)
        with open(out_dir+contig+"_msmc_input.log", 'w') as log:
            process = subprocess.Popen(cmd4, shell=True, stdout=log)
            # Set CPU affinity
            process_pid = process.pid
            os.sched_setaffinity(process_pid, cpu_affinity)
    process.wait()
    # Run MSMC2 combining information from all contigs
    cmd5 = " ".join(["msmc2_Linux", out_dir+"*_msmc_input.txt -o", out_dir+popid+"_msmc2"])
    print(cmd5)
    with open(popid+"_msmc2.log", 'w') as log:
        p=subprocess.Popen(cmd5,stdout=log, shell=True)
def psmc(ref_genome, contigs, popid, pop_ind, vcf, out_dir, mu, gen_time, num_cpus=4):
    # Get the available CPUs
    available_cpus = list(range(multiprocessing.cpu_count()))
    # Choose the first 'num_cpus' CPUs
    cpu_affinity = available_cpus[:num_cpus]

    #TODO
    # for sample in pop_ind:
    #     cmd1 = " ".join(["bcftools consensus -I", vcf, "-s", sample, "-f", ref_genome, "-o", out_dir+"consensus_"+sample+".fa"])
    #     with open(out_dir+sample+"_consensus.log", 'w') as log:
    #         process = subprocess.Popen(cmd1, shell=True, stdout=log)
    #         # Set CPU affinity
    #         process_pid = process.pid
    #         os.sched_setaffinity(process_pid, cpu_affinity)
    # process.wait()
    for sample in pop_ind:
        cmd2 = " ".join(["bcftools consensus -I", vcf, "-s", sample, "-f", ref_genome, "|",
        # read from stdin
        "seqtk seq -F '#'", "-", "|",
        "bgzip >", out_dir+"consensus_"+sample+".fq.gz", ";",
        "fq2psmcfa -q1", out_dir+"consensus_"+sample+".fq.gz", ">", out_dir+sample+"_diploid.psmcfa", ";"
        "psmc -N30 -t15 -r5 -p '10+22*2+4+6' -o", out_dir+sample+".psmc", out_dir+sample+"_diploid.psmcfa"])
        print(cmd2)
        with open(out_dir+sample+"_consensus.log", 'w') as log:
            process = subprocess.Popen(cmd2, shell=True, stdout=log)
            # Set CPU affinity
            process_pid = process.pid
            os.sched_setaffinity(process_pid, cpu_affinity)
    process.wait()
    cmd3 = " ".join(["cat", out_dir+"*.psmc >", out_dir+"combined.psmc.final", ";"
    "psmc_plot.pl", "-g", gen_time, "-x", "10", "-u", mu, "-M '"+",".join(pop_ind)+"'", popid, out_dir+"combined.psmc.final" ,";",
    "mv", popid+".*.par", out_dir, "; mv", popid+"*.eps", out_dir])
    print(cmd3)
    with open(out_dir+popid+"_psmc_combine.log", 'w') as log:
        p=subprocess.Popen(cmd3,stdout=log, shell=True)
    p.wait()

def smcpp(contigs, popid, pop_ind, vcf, out_dir, mu, gen_time):
    POP = popid+":"+",".join(pop_ind)
    if len(contigs) == 0:
        print("Error! No contigs to use! Make sure the threshold matches your data.")
        with open(out_dir+popid+"_model.final.json", 'w') as output:
            output.write("There was an error with SMC++. Please check the logs.")
        exit(0)

    for contig in contigs:
        cmd1 = ["smc++", "vcf2smc", vcf, out_dir+popid+"_"+contig+".smc.gz", contig, POP]
        print("LOG", out_dir+popid+"_"+contig+"_vcf2smc.log")
        with open(out_dir+popid+"_"+contig+"_vcf2smc.log", 'w') as log:
            p=subprocess.Popen(cmd1,stdout=log)
    p.wait()
    cmd2 = "".join(["smc++ estimate -o ", out_dir, popid, "_model.final.json ", mu, " ", out_dir, popid, "_*.smc.gz"])
    output2 = subprocess.check_output([cmd2], shell=True)
    output2 = output2.decode(encoding="utf-8")
    print(output2)
    cmd3 = "".join(["smc++ plot ", out_dir, popid, "_inference.png ", out_dir, popid, "_model.final.json/model.final.json -g", str(gen_time)," -c "])
    os.system(cmd3)

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Creates the input files to run the inferences softwares and run the softwares
"""

import re
import os
import subprocess
import itertools

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
def dadi_output_parse(dadi_output_file):
    all_vals = []
    with open(dadi_output_file) as dadi_output:
        for line in dadi_output:
            #print(line)
            line_parsed = line.strip().split()
            #print(line_parsed)
            ite = int(line_parsed[0])
            logL = float(line_parsed[2][:-1])
            all_vals.append([ite, logL, eval(line.strip().split("array")[-1][1:-1])])
    return all_vals




def dadi_inf(popid,out_dir,out_dir_d,dict,p,lower_bound,upper_bound,p0,mu,L,gen):
    
    import dadi
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


    

    





def input_stairwayplot2(popid, nseq, L, whether_folded, SFS, mu, year_per_generation, stairway_plot_dir, output_path, temp_blueprint):
    #writes the input file to run stairwayplot2
    #popid: name of the population
    #nseq: number of sequences (2*n for diploids)
    #L: total number of sites?
    #whether_folded: if the vcf is phased or not
    #SFS: sfs as a list without the monomorphic sites
    #mu: mutation rate
    #year_per_generation: generation time
    #stairway_plot_dir: path to stairwayplot2 (ce sera pas utile si on fait un conda j'imagine)
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
    cmd1 = "".join(["qsub -cwd -V -N stairwayplot2_", popid, " -o /home/tforest/work/birdsdemography/out.err/stairwayplot2_", popid, ".out -e /home/tforest/work/birdsdemography/out.err/stairwayplot2_", \
                   popid, ".err -q short.q -b y 'java -cp ", path_to_stairwayplot2, " Stairbuilder ", out_dir, popid, ".blueprint'"])
    #print(cmd1)
    output = subprocess.check_output([cmd1], shell=True) 
    job_submit_output = output.decode(encoding="utf-8")
    submitted_job_id = job_submit_output.split(" ")[2] #liste de l'ID des jobs
    #print(submitted_job_id)
    cmd2 = "".join(["qsub -cwd -V -N stairwayplot2_inf_", popid, " -o /home/tforest/work/birdsdemography/out.err/stairwayplot2_inf_", popid, ".out -e /home/tforest/work/birdsdemography/out.err/stairwayplot2_inf_", \
                   popid, ".err -q short.q -hold_jid ", submitted_job_id, " -b y 'bash ", out_dir, popid, ".blueprint.sh'"])
    #print(cmd2)
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


def smcpp(contigs, popid, pop_ind, vcf, out_dir, mu, gen_time):
    POP = popid+":"+",".join(pop_ind)
    #print(contigs)
    #subprocess.run(['/home/tforest/work/birdsdemography/smcp.sh','contigs','POP'],shell=True)
    jobs_list = []
    for contig in contigs:
        cmd1 = "".join(["qsub -cwd -V -N smcpp_", popid, "_", contig, " -o /home/tforest/work/birdsdemography/out.err/smcpp_", popid, "_", contig, \
                       ".out -e /home/tforest/work/birdsdemography/out.err/smcpp_", popid,"_", contig, \
                       ".err -q short.q -b y '/home/tforest/work/.conda/envs/smcpp/bin/smc++ vcf2smc ", vcf, " ", out_dir, popid, "_", contig, ".smc.gz ", contig, " ", POP, "'"])
        #print(cmd1)
        output = subprocess.check_output([cmd1], shell=True)
        job_submit_output = output.decode(encoding="utf-8")
        submitted_job_id = job_submit_output.split(" ")[2] #liste de l'ID des jobs
        #submitted_job_id = '143843'
        jobs_list.append(str(submitted_job_id))
    cmd2 = "".join(["qsub -cwd -V -N smcpp_inf_", popid, " -o /home/tforest/work/birdsdemography/out.err/smcpp_inf_", popid, \
                    ".out -e /home/tforest/work/birdsdemography/out.err/smcpp_inf_", popid, \
                    ".err -q short.q -hold_jid ", ",".join(jobs_list), " -b y '/home/tforest/work/.conda/envs/smcpp/bin/smc++ estimate -o ", out_dir, popid, "_model.final.json ", mu, \
                    " ", out_dir, popid, "_*.smc.gz'"])
    #print(cmd2)
    output2 = subprocess.check_output([cmd2], shell=True)
    job_submit_output2 = output2.decode(encoding="utf-8")
    submitted_job_id2 = job_submit_output2.split(" ")[2] #liste de l'ID des jobs
    cmd3 = "".join(["qsub -cwd -V -N smcpp_plot_", popid, " -o /home/tforest/work/birdsdemography/out.err/smcpp_plot_", popid, \
                ".out -e /home/tforest/work/birdsdemography/out.err/smcpp_plot_", popid, \
                ".err -q short.q -hold_jid ", submitted_job_id2, " -b y '/home/tforest/work/.conda/envs/smcpp/bin/smc++ plot ", \
                out_dir, popid, "_inference.png ", out_dir, popid, "_model.final.json/model.final.json -g", str(gen_time)," -c '"])
    os.system(cmd3)
    #print(cmd3)























# -*- coding: utf-8 -*-
"""
usage: project.py [-h] [--vcf VCF] [--n N] [--folded] [--transformed] [--sfs]
                  [--plot] [--stairwayplot2]
                  [--path_to_stairwayplot2 PATH_TO_STAIRWAYPLOT2]
                  [--plot_swplot2] [--summary_file SUMMARY_FILE]
                  [--name_plot NAME_PLOT] [--dadi] [--mu MU]
                  [--gen_time GEN_TIME]
                  popid

Computes the sfs from the vcf and runs demography inference softwares.

positional arguments:
  popid                 name of the population

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF             path to the vcf
  --n N                 number of individuals
  --folded              if the vcf is phased
  --transformed         if the sfs has to be transformed
  --sfs                 to compute the sfs
  --plot                to plot the sfs
  --stairwayplot2       to run stairwayplot2: the sfs must not be transformed,
                        the path to stairwayplot2 (--path-to-stairwayplot2),
                        the mutation rate (--mu) and the generation time
                        (--gen_time) must be specified
  --path_to_stairwayplot2 PATH_TO_STAIRWAYPLOT2
                        path to stairwayplot2
  --plot_swplot2        plot sariwayplot2 (by default, output file names
                        stairwayplot_plot_popid, or can be specifed with
                        --name_plot
  --summary_file SUMMARY_FILE
                        path to the plot summary file
  --name_plot NAME_PLOT
                        optional name of the sariwayplot2 plot
  --dadi                to run dadi : the sfs must not be transformed
  --mu MU               mutation rate (per site per generation)
  --gen_time GEN_TIME   generation time (in years)
"""


## MODULES 

#import os
import gzip
import itertools 
import matplotlib.pyplot as plt
import argparse
import numpy as np


## FONCTIONS

def compute_sfs(vcf, n, folded):
    #to compute the sfs from the vcf
    #vcf: path to the vcf file (str)
    #n: number of individuals (int)
    #folded: if the vcf is phased or not (bool)
    nsites = 0
    if folded: 
        sfs = [0] * n
    else:
        sfs = sfs = [0] * (2 * n - 1) 
    with gzip.open(vcf,  mode='rt') as file: 
        line = file.readline()
        while line != "": 
            if line[0] != "#" and "/." not in line and "./" not in line:   #il faudrait mettre une expression régulière?
                nsites += 1
                split_line = line.split("\t") 
                if "," not in split_line[4]: #we only keep the bi-allelique sites
                    gen = list(itertools.chain.from_iterable([[int(i[0]), int(i[2])] for i in split_line[-(n):]])) #we get the genotypes
                    #gen a list of 1 and 0
                    if folded:
                        count=min(sum(gen), 2*n-sum(gen)) #if folded, we count the minor allele
                    else:
                        count=sum(gen) #if not folded, we count the alternate allele (1)
                    if count != 0 and (folded) or (not folded and count != 2*n): #we remove the monomorphic sites
                        #print(gen)
                        sfs[count - 1] = sfs[count - 1] + 1 #the sfs is incremented
            line = file.readline()
    return(sfs, nsites)


def transform_sfs(sfs, n, folded):
    #to transform the vcf
    #formula from Elise thesis
    if folded:
        return [round((i + 1) * (2 * n - ( i + 1)) * sfs[i] / (2 * n)) for i in range(len(sfs))]
    else:
        return [round((i + 1) * sfs[i]) for i in range(len(sfs))]


def plot_sfs(sfs, n, folded, transformed, popid):
    #to output the plot
    if folded:
        x = np.arange(1, n+1)
    else:
        x = np.arange(1, 2*n)
    if transformed:
        plot_title = 'Site Frequency Spectrum - transformed'
        name_output = 'sfs_'+str(popid)+'_transformed.png'
    else:
        plot_title = 'Site Frequency Spectrum'
        name_output = 'sfs_'+str(popid)+'.png'
    plt.bar(x, sfs)
    plt.title(plot_title)
    plt.xlabel('# of alternate alleles')
    plt.ylabel('# sites')
    plt.savefig(name_output)



def input_stairwayplot2(popid, nseq, L, whether_folded, SFS, mu, year_per_generation, stairway_plot_dir):
    #writes the input file to run stairwayplot2
    #popid: name of the population
    #nseq: number of sequences (2*n for diploids)
    #L: total number of sites?
    #whether_folded: if the vcf is phased or not
    #SFS: sfs as a list without the monomorphic sites
    #mu: mutation rate
    #year_per_generation: generation time
    #stairway_plot_dir: path to stairwayplot2 (ce sera pas utile si on fait un conda j'imagine)
    locals()['project_dir'] = popid #where the output of stairwayplot2 will be
    locals()['plot_title'] = popid #name of the plots output by stairwayplot2
    with open("template.blueprint", "r") as temp, open(str(popid)+".blueprint","w") as out_file:
        line = temp.readline()
        while line != '':
            if line.split(':')[0] in locals().keys():
                if line.split(':')[0] == 'SFS':
                    out_file.write('sfs: ' + ' '.join(map(str,sfs)) + '\n') 
                else:
                    out_file.write(line.split(':')[0]+': '+ str(locals().get(line.split(':')[0]))+'\n')
            elif line.split(':')[0] == "nrand" :
                out_file.write('nrand: '+str(int((nseq-2)/4))+' '+ str(int((nseq-2)/2))+' '+ str(int((nseq-2)*3/4))+' '+ str(int(nseq-2))+'\n')
            else:
                out_file.write(line)
            line = temp.readline()
            
            
            
def plot_stairwayplot2(popid, summary_file, output_name):
    Ne_med=[]
    Ne1=[]
    Ne3=[]
    T=[]
    with open(summary_file) as input_file:
        line = input_file.readline()
        line = input_file.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T: 
            Ne_med.append(float(line.split('\t')[6]))
            Ne1.append(float(line.split('\t')[7]))
            Ne3.append(float(line.split('\t')[8]))
            T.append(float(line.split('\t')[5]))
            line = input_file.readline()
    plt.plot(T,Ne_med,color="blue")
    plt.plot(T,Ne1,color="grey")
    plt.plot(T,Ne3,color="grey")
    plt.title(popid)
    plt.xlabel('Time (in years)')
    plt.ylabel('Ne')
    plt.savefig(output_name)
    


def input_dadi (popid, sfs, folded, n, mask = True):
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
    with open(popid+'_dadi_input.txt', 'w') as f:
        f.write(dim+'\n'+sfs_write+"\n"+["0","1"][mask])


    

## SCRIPT

parser = argparse.ArgumentParser(description='Computes the sfs from the vcf and runs demography inference softwares.')

#mandatory arguments
parser.add_argument("popid", help="name of the population")

#optional arguments
parser.add_argument("--vcf", help="path to the vcf")
parser.add_argument("--n", type=int, help="number of individuals")
parser.add_argument("--folded", help="if the vcf is phased", action = "store_true")
parser.add_argument("--transformed", help = "if the sfs has to be transformed", action = "store_true")
parser.add_argument("--sfs", help = "to compute the sfs", action = "store_true")
parser.add_argument("--plot", help = "to plot the sfs", action = "store_true")
parser.add_argument("--stairwayplot2", help = "to run stairwayplot2: the sfs must not be transformed, the path to stairwayplot2 (--path-to-stairwayplot2), the mutation rate (--mu) and the generation time (--gen_time) must be specified", action = "store_true")
parser.add_argument("--path_to_stairwayplot2", help = "path to stairwayplot2")
parser.add_argument("--plot_swplot2", help = "plot sariwayplot2 (by default, output file names stairwayplot_plot_popid, or can be specifed with --name_plot ", action = "store_true")
parser.add_argument("--summary_file", help = "path to the plot summary file")
parser.add_argument("--name_plot", help = "optional name of the sariwayplot2 plot")
parser.add_argument("--dadi", help = "to run dadi : the sfs must not be transformed", action = "store_true")
parser.add_argument("--mu", type=float, help="mutation rate (per site per generation)")
parser.add_argument("--gen_time", type=int, help="generation time (in years)")

args = parser.parse_args()


#vcf = "/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/data/Historical_Merged.vcf.gz"
#n = 5
#folded = True
#transformed = False
#popid = 'hiron_contemp'
#mu = 10e-8
#gen_time = 2
#plot = True
#stairwayplot2 = True
#path_to_stairwayplot2 = "/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/scripts/stairway_plot_v2.1.1/stairway_plot_es"
#dadi = True
#nrand=7
#plot_swplot2 = True
#summary_file = "/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/scripts/hiron_contemp/hiron_contemp.final.summary"


if args.sfs:
    sfs = compute_sfs(args.vcf, args.n, args.folded)[0]
    nsites = compute_sfs(args.vcf, args.n, args.folded)[1]
    if args.transformed:
        sfs = transform_sfs(sfs, args.n, args.folded)
        with open(args.popid+"_sfs_transformed.txt", 'w') as sfs_out:
            for i in sfs:
                sfs_out.write(str(i)+" ")
    else:
        with open(args.popid+"_sfs.txt", 'w') as sfs_out:
            for i in sfs:
                sfs_out.write(str(i)+" ")

if args.plot:
    plot_sfs(sfs, args.n, args.folded, args.transformed, args.popid)
    
if args.stairwayplot2:
    if args.sfs:
        input_stairwayplot2(popid = args.popid, nseq = 2*args.n, L= nsites, whether_folded = args.folded, 
                            SFS = sfs, mu = args.mu, year_per_generation = args.gen_time, stairway_plot_dir=args.path_to_stairwayplot2)
    else:
        sfs = compute_sfs(args.vcf, args.n, args.folded)[0]
        nsites = compute_sfs(args.vcf, args.n, args.folded)[1]
        input_stairwayplot2(popid = args.popid, nseq = 2*args.n, L= nsites, whether_folded = args.folded, 
                            SFS = sfs, mu = args.mu, year_per_generation = args.gen_time, stairway_plot_dir=args.path_to_stairwayplot2)

if args.plot_swplot2:
    if args.name_plot:
        plot_stairwayplot2(args.popid, args.summary_file, args.name_plot)
    else:
        plot_stairwayplot2(args.popid, args.summary_file, "stairwayplot2_plot_"+args.popid+".png")

    
if args.dadi:
    if args.sfs:
        input_dadi(args.popid, sfs, args.folded, args.n, True)
    else:
        sfs = compute_sfs(args.vcf, args.n, args.folded)[0]
        input_dadi(args.popid, sfs, args.folded, args.n, True)


            













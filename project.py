# -*- coding: utf-8 -*-
"""
usage: project.py [-h] [--transformed] [--plot] [--stairwayplot2]
                  [--path_to_stairwayplot2 PATH_TO_STAIRWAYPLOT2] [--dadi]
                  [--mu MU] [--gen_time GEN_TIME]
                  vcf n popid folded

Computes the sfs from the vcf and runs demography inference softwares. By
default, outputs the non-transformed sfs in a text file named {popid}_sfs.txt

positional arguments:
  vcf                   path to the vcf
  n                     number of individuals in the pop
  popid                 name of the population
  folded                either the vcf is phased or not

optional arguments:
  -h, --help            show this help message and exit
  --transformed         if the sfs has to be transformed
  --plot                to plot the sfs
  --stairwayplot2       to run stairwayplot2: the sfs must not be transformed,
                        the path to stairwayplot2 (--path-to-stairwayplot2),
                        the mutation rate (--mu) and the generation time
                        (--gen_time) must be specified
  --path_to_stairwayplot2 PATH_TO_STAIRWAYPLOT2
                        path to stairwayplot2
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

parser = argparse.ArgumentParser(description='Computes the sfs from the vcf and runs demography inference softwares. By default, outputs the non-transformed sfs in a text file named {popid}_sfs.txt')

#mandatory arguments
parser.add_argument("vcf", help="path to the vcf")
parser.add_argument("n", type=int, help="number of individuals in the pop")
parser.add_argument("popid", help="name of the population")
parser.add_argument("folded", type=bool, help="either the vcf is phased or not")

#optional arguments
parser.add_argument("--transformed", help = "if the sfs has to be transformed", action = "store_true")
parser.add_argument("--plot", help = "to plot the sfs", action = "store_true")
parser.add_argument("--stairwayplot2", help = "to run stairwayplot2: the sfs must not be transformed, the path to stairwayplot2 (--path-to-stairwayplot2), the mutation rate (--mu) and the generation time (--gen_time) must be specified", action = "store_true")
parser.add_argument("--path_to_stairwayplot2", help = "path to stairwayplot2")
parser.add_argument("--dadi", help = "to run dadi : the sfs must not be transformed", action = "store_true")
parser.add_argument("--mu", type=float, help="mutation rate (per site per generation)")
parser.add_argument("--gen_time", type=int, help="generation time (in years)")

args = parser.parse_args()


#vcf = "/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/data/Historical_Merged.vcf.gz"
#n = 5
#folded = True
#transformed = False
#popid = 'hirundo_historical'
#mu = 10e-8
#gen_time = 2
#plot = True
#stairwayplot2 = True
#path_to_stairwayplot2 = "/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/scripts/stairway_plot_v2.1.1/stairway_plot_es"
#dadi = True


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
    input_stairwayplot2(popid = args.popid, nseq = 2*args.n, L= nsites, whether_folded = args.folded, 
                            SFS = sfs, mu = args.mu, year_per_generation = args.gen_time, stairway_plot_dir=args.path_to_stairwayplot2)
    
if args.dadi:
    input_dadi(args.popid, sfs, args.folded, args.n, True)


            













# -*- coding: utf-8 -*-
"""
Script pour lire vcf et calculer le sfs.
Ligne de commande : python sfs.py file folded standardized
où file est le nom du fichier vcf.gz, folded est un booleen indiquant si sfs plié 
ou non, et standardized est un booléen pour avoir le sfs standardisé ou non.
"""

#Importation des modules
import os
import gzip
import itertools 
import matplotlib.pyplot as plt
import sys
import numpy as np



## FONCTIONS

def sfs(vcf, n, folded):
    #on initialise un sfs
    if folded: 
        sfs = [0] * n
    else:
        sfs = sfs = [0] * (2 * n - 1) 
    with gzip.open(vcf,  mode='rt') as file: 
        line = file.readline()
        while line != "": #tant qu'il reste des lignes a parcourir
            if line[0] != "#" and "/." not in line and "./" not in line:   #il faudrait mettre une expression régulière
                split_line = line.split("\t") 
                if "," not in split_line[4]: #on garde que les sites bi-alléliques (à ajouter dans express. reg. ?)
                    gen = list(itertools.chain.from_iterable([[int(i[0]), int(i[2])] for i in split_line[-(n):]])) #on recupere les genotypes 
                    #gen est une liste d'entier O ou 1
                    if folded:
                        count=min(sum(gen), 2*n-sum(gen)) #on compte l'allele minoritaire
                    else:
                        count=sum(gen) # on compte l'allele alternatif
                    if count != 0 and (folded) or (not folded and count != 2*n): #on retire les sites monomorphiques
                        #print(gen)
                        sfs[count - 1] = sfs[count - 1] + 1 #on incrémente 
            line = file.readline()
    return(sfs)


def transform_sfs(sfs, n, folded):
    if folded:
        return [(i + 1) * (2 * n - ( i + 1)) * sfs[i] / (2 * n) for i in range(len(sfs))]
    else:
        return [(i + 1) * sfs[i] for i in range(len(sfs))]


def plot_sfs(sfs, n, folded, transformed, popid):
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



def blueprint_stairwayplot2(popid, nseq, L, whether_folded, SFS, mu, year_per_generation):
    locals()['project_dir'] = popid
    locals()['plot_title'] = popid
    print(locals().keys())
    with open("template.blueprint", "r") as temp, open(str(popid)+".blueprint","w") as out_file:
        line = temp.readline()
        while line != '':
            print(line.split(':')[0])
            if line.split(':')[0] in locals().keys():
                print("True")
                out_file.write(line.split(':')[0]+': '+ str(locals().get(line.split(':')[0]))+'\n')
            elif line.split(':')[0] == "nrand" :
                out_file.write('nrand :'+str(int((nseq-2)/4))+' '+ str(int((nseq-2)/2))+' '+ str(int((nseq-2)*3/4))+' '+ str(int(nseq-2))+'\n')
            else:
                out_file.write(line)
            line = temp.readline()



##############################################################################


os.chdir("/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/scripts")

## Arguments
#vcf, folded, transformed, n = sys.argv[1:5]


vcf = "/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/data/Contemporary_Merged_test.vcf.gz"
n = 6
folded = True
transformed = False
popid = 'hiron'
L = 100000
mu = 10e-8
gen_time = 2
plot = True
stairwayplot2 = True


sfs = sfs(vcf, n, folded)

if transformed:
    sfs = transform_sfs(sfs, n, folded)

if plot:
    plot_sfs(sfs, n, folded, transformed, popid)

if stairwayplot2:
    blueprint_stairwayplot2(popid = popid, nseq = 2*n, L= L, whether_folded = folded, 
                            SFS = sfs, mu = mu, year_per_generation = gen_time)














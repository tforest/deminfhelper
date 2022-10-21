# -*- coding: utf-8 -*-
"""
Script pour lire vcf et calculer le sfs.
Ligne de commande : python3 sfs.py file folded standardized
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

def compute_sfs(vcf, n, folded):
    #on initialise un sfs
    nsites = 0
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
                    nsites += 1
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
    return(sfs, nsites)


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



def input_stairwayplot2(popid, nseq, L, whether_folded, SFS, mu, year_per_generation):
    locals()['project_dir'] = popid
    locals()['plot_title'] = popid
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



def input_dadi (sfs, folded, n, mask = True):
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
    sfs_write=(''.join(map(str,sfs_write))) #concaténation de la liste un une chaîne de caractères
    sfs_write=sfs_write[0:len(sfs_write)-1]#retrait du dernier espace de la ligne
    #écriture du fichier :
    with open('dadi_input.txt', 'w') as f:
        f.write(dim+'\n'+sfs_write+"\n"+["0","1"][mask])



    

##############################################################################


## Arguments
#vcf, folded, transformed, n = sys.argv[1:5]

os.chdir("/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/scripts")

vcf = "/media/camille/Donnees/Documents/Etudes/ENS/M2/Projet/data/Historical_Merged.vcf.gz"
n = 5
folded = True
transformed = True
popid = 'hirundo_historical'
mu = 10e-8
gen_time = 2
plot = True
stairwayplot2 = True
dadi = True

sfs = compute_sfs(vcf, n, folded)[0]
nsites = compute_sfs(vcf, n, folded)[1]

if transformed:
    sfs = transform_sfs(sfs, n, folded)

if plot:
    plot_sfs(sfs, n, folded, transformed, popid)
    
if stairwayplot2:
    input_stairwayplot2(popid = popid, nseq = 2*n, L= nsites, whether_folded = folded, 
                            SFS = sfs, mu = mu, year_per_generation = gen_time)
    
if dadi:
    input_dadi(sfs, folded, n, True)
            













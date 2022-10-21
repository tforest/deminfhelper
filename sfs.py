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


<<<<<<< HEAD
## Arguments
name_file, folded, standardized = sys.argv[1:4]
#name_file, folded, standardized = "Contemporary_Merged_test.vcf", True, True
print(name_file, folded, standardized)


with gzip.open(name_file,  mode='rt') as file:
    line = file.readline()
    while line != "": #tant qu'il reste des lignes a parcourir
        if line[0:6] == "#CHROM": #on recupere le nombre de samples
            n = len(line.split("\t")[9:])
            sfs = [0] * (2 * n - 1) #on initialise un sfs
        if line[0] != "#" and "/." not in line and "./" not in line:  #on parcours le fichier jusqu'aux snps
        #il faudrait mettre une expression régulière
            split_line = line.split("\t")
            if "," not in split_line[4]: #on garde que les sites bi-alléliques (à ajouter dans express. reg. ?)
                gen = list(itertools.chain.from_iterable([[int(i[0]), int(i[2])] for i in split_line[-(n):]])) #on recupere les genotypes
                #gen est une liste d'entier O ou 1
                #sum(gen) est le nombre d'alleles alternatifs
                if sum(gen) != 2 * n and sum(gen) != 0: #on retires les sites monomorphiques
                    #print(gen)
                    sfs[sum(gen) - 1] = sfs[sum(gen) - 1] + 1 #on incrémente
        line = file.readline()

#on plie et/ou standardise si demande
if folded:
    sfs2 = [sfs[i-1]+sfs[-i] for i in range(1,n)]
    if standardized:
        sfs_std = [(i + 1) * (2 * n - ( i + 1)) * sfs2[i] / (2 * n) for i in range(len(sfs2))]
else:
    sfs2 = sfs
    if standardized:
        sfs_std = [(i + 1) * sfs2[i] for i in range(len(sfs2))]


#Outputs (only plots for now)
plt.bar(np.arange(1, len(sfs2) + 1), sfs2)
plt.title('Site Frequency Spectrum')
plt.xlabel('# of alternate alleles')
plt.ylabel('# sites')
plt.savefig('SFS_nonstand_'+name_file[:-7]+'.png')
plt.clf()

plt.bar(np.arange(1, len(sfs_std) + 1), sfs_std, )
plt.title('Site Frequency Spectrum - standardized')
plt.xlabel('# of alternate alleles')
plt.ylabel('# sites')
plt.savefig('SFS_stand_'+name_file[:-7]+'.png')

#Outputs: we'll have to take in arguments the software used next (dadi or stairwayplot) to format the sfs the correct way
=======

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
>>>>>>> 15b42af994a9b5f9108638f15b11ba0df2b016b8

###

#os.chdir("/BIRDS/VCF_HIR")


import os
import gzip
import itertools
import matplotlib.pyplot as plt
import sys
import numpy as np


def dadi (sfs,folded,n,mask):
    #N_pops est l'effectif d'une
    #sfs est une np.array avec chaque colonne correspondant à un effectif de variant et chaque ligne
    #correspondant à une population
    #mask = True si masked, False sinon
    #fold = True si folded, False sinon
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



dadi (sfs,folded,n,True)

[]

####

<<<<<<< HEAD
=======
#
# sfs_std=np.array([[0, 0, 0, 0, 0, 1, 0, 0, 1, 4, 2],[0, 0, 0, 0, 0, 2 ,5 , 2 ,3 , 4, 2]])
#
# Npops = [2*n,6]

# def dadi (sfs, folded, N_pops, mask):
#     #N_pops est une liste contenant l'effectif de chaque population
#     #sfs est une np.array avec chaque colonne correspondant à un effectif de variant et chaque ligne
#     #correspondant à une population
#     #mask = True si masked, False sinon
#     #fold = True si folded, False sinon
#     N_pops=[i+1 for i in N_pops]
#     #Première ligne :
#     dim=""
#     for i in N_pops:
#         if len(dim)==0: #Pour ne pas ajouter d'espace au début de la ligne
#             dim=dim+str(i)
#         else :
#             dim=dim+" "+str(i)#ajoute un espace et la dimension à la première ligne
#     dim =dim+" "+["unfolded","folded"][fold]
#     #ligne de sfs :
#     sfs_write=list(itertools.chain.from_iterable([["fs",list(sfs[:,i])," "] for i in range (np.shape(sfs)[1])]))
#     sfs_write=(''.join(map(str,sfs_write))) #concaténation de la liste un une chaîne de caractères
#     sfs_write=sfs_write[0:len(sfs_write)-1]#retrait du dernier espace de la ligne
#     #écriture du fichier :
#     with open('dadi_input.txt', 'w') as f:
#         f.write(dim)
#         f.write('\n')
#         f.write(sfs_write)
#         f.write("\n"+["0","1"][mask])
#
#
#
# dadi (sfs_std,folded,Npops,True)
#
>>>>>>> 15b42af994a9b5f9108638f15b11ba0df2b016b8




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














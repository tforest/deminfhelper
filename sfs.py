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
os.chdir("/BIRDS/VCF_HIR")

## Arguments
name_file, folded, standardized = sys.argv[1:4]
print(name_file, folded, standardized)


with open(name_file,  mode='rt') as file:
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

#
# sfs_std=np.array([[0, 0, 0, 0, 0, 1, 0, 0, 1, 4, 2],[0, 0, 0, 0, 0, 2 ,5 , 2 ,3 , 4, 2]])
#
# Npops = [2*n,6]
# def dadi (sfs,fold,N_pops,mask):
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

















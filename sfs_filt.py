#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
build_sfs (n, folded,sfs_ini, line='', sfs=[]) : compute the sfs
transform_sfs(sfs, n, folded) : transforms the sfs
"""

import itertools
import numpy as np

def build_sfs(n, folded, sfs_ini, line=[], sfs=[], pos_ind = None):
    if sfs_ini:
        if folded:
            return [0] * n
        else:
            return [0] * (2 * n - 1)
    else:
        geno_ind = [line[i] for i in pos_ind]
        gen = list(itertools.chain.from_iterable([[int(i[0]), int(i[2])] for i in geno_ind])) #we get the genotypes
        #gen a list of 1 and 0
        if folded:
            count=min(sum(gen), 2*n-sum(gen)) #if folded, we count the minor allele
        else:
            count=sum(gen) #if not folded, we count the alternate allele (1)
        if count != 0 and (folded) or (not folded and count != 2*n): #we remove the monomorphic sites
            #print(gen)
            sfs[count - 1] = sfs[count - 1] + 1 #the sfs is incremented
        return(sfs)


def transform_sfs(sfs, n, folded):
    #to transform the vcf
    #formula from Elise thesis
    if folded:
        return [round((i + 1) * (2 * n - ( i + 1)) * sfs[i] / (2 * n)) for i in range(len(sfs))]
    else:
        return [round((i + 1) * sfs[i]) for i in range(len(sfs))]


def distrib_GQ(GQ_pop, line = [], pos_ind = None): #PL is a dict
    samples = [line[i] for i in pos_ind]
    gq_line = [i[4:] for i in samples] #we get the genotypes
    gq_line = [i.split(",") for i in gq_line]
    for sublist in gq_line:
        gq2 = [int(i) for i in sublist]
        gq2 = [i - np.min(gq2) for i in gq2]
        gq2.remove(0)
        min = np.min(gq2)
        bin = min - int(str(min)[-1])
        if bin in GQ_pop.keys():
            GQ_pop[bin] = GQ_pop[bin]+1
        else:
            GQ_pop[bin] = 1
    return(GQ_pop)



def check_GQ(n_ind, thresh, line = []): #PL is a dict
    gq_line = [i[4:] for i in line[9:n_ind+10]] #we get the genotypes
    gq_line = [i.split(",") for i in gq_line]
    for sublist in gq_line:
        gq2 = [int(i) for i in sublist]
        gq2 = [i - np.min(gq2) for i in gq2]
        gq2.remove(0)
        min = np.min(gq2)
        if min < thresh:
            return(False)
    return(True)







# -*- coding: utf-8 -*-
"""
parsing() : parses the vcd, compute the sfs if SFS=True, output the msmc input file
if MSMC=False
"""

import gzip
from inferences import *
from sfs_filt import *
import re


def parsing(PARAM, SFS = False, GQ = False, SMCPP = False):      
    with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
        SFS_dict = {}
        GQ_dict = {}
        contigs = []
        if SFS:
            # we initialize a sfs for each population
            for p in PARAM["name_pop"]:
                SFS_dict[p] = build_sfs(n=PARAM["n_"+str(p)], folded=PARAM["folded"], sfs_ini=True)
        if GQ:
            for p in PARAM["name_pop"]:
                        GQ_dict[p] = {}
        if SMCPP:
            contigs = []
        line = vcf.readline()
        # we read all the lines of the vcf
        while line != "":
            if line[0:8] == "##contig":
                if int(re.split('[=,]', line)[-1][:-2]) >= 100000:
                    contigs.append(re.split('[=,]', line)[2])
            if line[0:6] == "#CHROM":
                for p in PARAM["name_pop"]:
                    pos_ind_pop = []
                    for ind in PARAM[p]:
                        pos_ind_pop.append(line[:-1].split("\t").index(ind))
                    PARAM["pos_"+p] = pos_ind_pop
            if line[0] != "#" and "/." not in line and "./" not in line and "," not in line.split("\t")[4]:    #we only keep the bi-allelique sites
                split_line = line.split("\t")
                if SFS:
                    for p in PARAM["name_pop"]:
                        SFS_dict[p] = build_sfs(n=PARAM["n_"+p], folded=PARAM["folded"],  sfs_ini = False, \
                                line = split_line, sfs = SFS_dict[p], pos_ind = PARAM["pos_"+p])
                if GQ:
                    for p in PARAM["name_pop"]:
                        GQ_dict[p] = distrib_GQ(GQ_pop = GQ_dict[p], line = split_line, pos_ind = PARAM["pos_"+p])
            line = vcf.readline()

    return SFS_dict, GQ_dict, contigs




def parsing_filt(PARAM, thresh = 0, SFS = False, GQ = False, SMCPP = False):      
    with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
        nsites = 0
        SFS_dict = {}
        GQ_dict = {}
        contigs = []
        if SFS:
            outfile = open(PARAM["out_dir"]+PARAM["out_dir"].split("/")[-2]+".vcf", "w")
            # we initialize a sfs for each population
            for p in PARAM["name_pop"]:
                SFS_dict[p] = build_sfs(n=PARAM["n_"+str(p)], folded=PARAM["folded"], sfs_ini=True)
        if GQ:
            for p in PARAM["name_pop"]:
                        GQ_dict[p] = {}
        if SMCPP:
            contigs = []
        line = vcf.readline()
        # we read all the lines of the vcf
        while line != "":
            if SFS:
                if line[0:1] == "#":
                    outfile.write(line)
            if line[0:8] == "##contig":
                if int(re.split('[=,]', line)[-1][:-2]) >= 100000:
                    contigs.append(re.split('[=,]', line)[2])
            if line[0:6] == "#CHROM":
                for p in PARAM["name_pop"]:
                    pos_ind_pop = []
                    for ind in PARAM[p]:
                        pos_ind_pop.append(line[:-1].split("\t").index(ind))
                    PARAM["pos_"+p] = pos_ind_pop
            if SFS or GQ:
                if line[0] != "#" and "/." not in line and "./" not in line and "," not in line.split("\t")[4]:    #we only keep the bi-allelique sites
                    split_line = line.split("\t")
                    tot_ind = 0
                    for p in PARAM["name_pop"]:
                        tot_ind = tot_ind + PARAM["n_"+p]
                    if check_GQ(n_ind = tot_ind, line = split_line, thresh = thresh):
                        nsites += 1
                        if SFS:
                            outfile.write(line)
                            for p in PARAM["name_pop"]:
                                SFS_dict[p] = build_sfs(n=PARAM["n_"+p], folded=PARAM["folded"],  sfs_ini = False, \
                                        line = split_line, sfs = SFS_dict[p], pos_ind = PARAM["pos_"+p])
                        if GQ:
                            for p in PARAM["name_pop"]:
                                GQ_dict[p] = distrib_GQ(GQ_pop = GQ_dict[p], line = split_line, pos_ind = PARAM["pos_"+p])
            line = vcf.readline()

    return SFS_dict, GQ_dict, contigs, nsites











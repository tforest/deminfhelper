#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

"""

## IMPORT MODULES

import gzip
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

## CMD ARGUMENTS 

parser = argparse.ArgumentParser(description='Computes the sfs from the vcf and runs demography inference softwares.')

#mandatory arguments
parser.add_argument("config_file", help="path to the configuration file")

#optional arguments
#SFS
parser.add_argument("--sfs", help = "to compute the sfs", action = "store_true")
parser.add_argument("--sfs_transformed", help = "to normalize the sfs", action = "store_true")
parser.add_argument("--plot_sfs", help = "to plot the sfs", action = "store_true")
#Stairwayplot2
parser.add_argument("--stairwayplot2", help = "to run stairwayplot2", action = "store_true")
parser.add_argument("--plot_stairwayplot2", help = "to run stairwayplot2", action = "store_true")
#Dadi
parser.add_argument("--dadi", help = "to run dadi: the sfs must not be transformed", action = "store_true")
#MSMC
parser.add_argument("--msmc", help = "to run msmc: the sfs must not be transformed", action = "store_true")
#PL distribution
parser.add_argument("--gq_distrib", help = "to compute the GQ (genotype quality) distribution", action = "store_true")
#SMCPP
parser.add_argument("--smcpp", help = "run smcpp", action = "store_true")
parser.add_argument("--plot_smcpp", help = "to plot smcpp inference", action = "store_true")
parser.add_argument("--Gplot", help = "to plot all inferences on the same graph", action = "store_true")


args = parser.parse_args()

## CONFIG FILE
param = {}
with open(args.config_file, "rt") as config:
    line=config.readline()
    while line != "":
        if line[0] != "#":
                param[line[:-1].split(": ")[0]] = line[:-1].split(": ")[1]
        line = config.readline()


param["folded"]=bool(param["folded"])
#param["transformed"]=bool(param["transformed"])
param["name_pop"] = param["name_pop"].split(",")
param["npop"]=int(param["npop"])
for p in param["name_pop"]:
    param[p] = param[p].split(",")
    param["n_"+p] = len(param[p])


#for key in param.keys():
#    print(key,": ",param[key])


## CREATING DIRECTORIES
if not os.path.exists(param["out_dir"]):
    os.makedirs(param["out_dir"])

if args.smcpp or args.stairwayplot2 or args.dadi:
    if not os.path.exists(param["final_out_dir"]):
        os.makedirs(param["final_out_dir"])


## IMPORT FUNCTIONS
from parsing import *
from inferences import *
from sfs import *
from plots import *


# Compute the SFS

if args.sfs or args.gq_distrib or args.smcpp:
    res_pars = parsing(PARAM = param, SFS = args.sfs, SMCPP = args.smcpp, GQ = args.gq_distrib)

if args.sfs:
    if not os.path.exists(param["out_dir_sfs"]):
        os.makedirs(param["out_dir_sfs"])
    SFS_dict = res_pars[0]
    for p in param["name_pop"]:
        with open(param["out_dir_sfs"]+"SFS_"+p+".txt", 'w') as sfs_out:
            sfs_out.write(",".join(map(str, SFS_dict[p]))+"\n")
        if args.plot_sfs:
            plot_sfs(sfs = SFS_dict[p], plot_title = "SFS "+p, output_file = param["out_dir_sfs"]+p+".png")


if args.sfs_transformed:
    if args.sfs == False:
        SFS_dict = {}
        for p in param["name_pop"]:
            if "path_to_sfs_"+p not in param.keys():
                print("--sfs flag or path_to_sfs missing")
            else:
                with open(param["path_to_sfs_"+p], "rt") as sfs:
                    line=sfs.readline()
                    while line != "":
                        SFS_dict[p] = [int(i) for i in line[:-1].split(",")]
                        line = sfs.readline()
    SFS_dict_trans = {}
    for p in param["name_pop"]:
        SFS_dict_trans[p] = transform_sfs(sfs = SFS_dict[p], n = param["n_"+p], \
                folded = param["folded"])
        with open(param["out_dir_sfs"]+"SFS_"+p+"_transformed.txt", 'w') as sfs_out:
            sfs_out.write(",".join(map(str, SFS_dict_trans[p])))
        if args.plot_sfs:
            plot_sfs(sfs = SFS_dict_trans[p], plot_title = "SFS (transformed) "+p, output_file = param["out_dir_sfs"]+p+"_trans.png")



# Run Stairwayplot2
if args.stairwayplot2:
    if not os.path.exists(param["out_dir_stairwayplot2"]):
        os.makedirs(param["out_dir_stairwayplot2"])
    if args.sfs == False:
        if "path_to_sfs_"+p not in param.keys():
            print("--sfs flag or path_to_sfs missing")
        else:
            SFS_dict = {}
#            with open(param["path_to_sfs_"+p], "rt") as sfs:
#                line=sfs.readline()
#                while line != "":
#                    SFS_dict[line[:-1].split(": ")[0]] = [int(i) for i in line[:-1].split(": ")[1].split(",")]
#                    line = sfs.readline()
            with open(param["path_to_sfs_"+p], "rt") as sfs:
                line=sfs.readline()
                while line != "":
                    SFS_dict[p] = [int(i) for i in line[:-1].split(",")]
                    line = sfs.readline() 

    print(SFS_dict)
    for p in param["name_pop"]:
        input_stairwayplot2(popid = p, nseq = param["n_"+p]*2, L= param["L"], whether_folded = param["folded"], \
                        SFS = SFS_dict[p] , mu = param["mut_rate"], year_per_generation = param["gen_time"], \
                        stairway_plot_dir = param["path_to_stairwayplot2"], output_path = param["out_dir_stairwayplot2"], \
                        temp_blueprint = param["blueprint_template"])
        run_stairwayplot2(popid = p, out_dir = param["out_dir_stairwayplot2"], path_to_stairwayplot2 = param["path_to_stairwayplot2"])
    if args.plot_stairwayplot2:
        for p in param["name_pop"]:
            plot_stairwayplot2(popid = p, summary_file = "".join([param["out_dir_stairwayplot2"], p, "/", p,".final.summary"]), \
                           out_dir = param["final_out_dir"])


if args.plot_stairwayplot2 and args.stairwayplot2==False:
    if not os.path.exists(param["out_dir_stairwayplot2"]):
        os.makedirs(param["out_dir_stairwayplot2"])
    for p in param["name_pop"]:
        if "".join(["summary_file_stw_",p]) not in param.keys():
            print("path to the population final summary file missing")
        else:
            plot_stairwayplot2(popid = p, summary_file = param["summary_file_stw_"+p], \
                           out_dir = param["final_out_dir"])




# Run dadi
if args.dadi:
    import dadi
    if not os.path.exists(param["out_dir_dadi"]):
        os.makedirs(param["out_dir_dadi"])
    if args.sfs == False:
        SFS_dict = {}
        for p in param["name_pop"]:
            if "path_to_sfs_"+p not in param.keys():
                print("--sfs flag or path_to_sfs missing")
            else:
                with open(param["path_to_sfs_"+p], "rt") as sfs:
                    line=sfs.readline()
                    while line != "":
                        SFS_dict[p] = [int(i) for i in line[:-1].split(",")]
                        line = sfs.readline() 
        
    for p in param["name_pop"]:
        print(SFS_dict)
        DICT = {i : j for i, j in zip([i for i in range(0,2*len(SFS_dict[p]))], SFS_dict[p]+[0]*len(SFS_dict[p]))}
        input_dadi(popid = p, sfs = SFS_dict[p], folded = param["folded"], n = param["n_"+p], \
                   mask = True, out_dir = param["out_dir_dadi"])
        print (type(param["lower_bound"]))
        dadi_inf(popid = p, out_dir_d = param["out_dir_dadi"],out_dir=param["out_dir"],dict=DICT,p=p,lower_bound=[float(x) for x in str.split(param["lower_bound"],',')],upper_bound=[float(x) for x in str.split(param["upper_bound"],',')],p0=[float(x) for x in str.split(param["p0"],',')],mu=param["mut_rate"],L=param["L"],gen=param["gen_time"])
        all_vals_recent = dadi_output_parse(param["out_dir_dadi"]+"output_"+p+".dadi")
        print_dadi_output_two_epochs(T_scaled_gen=param["out_dir_dadi"]+"popt_"+p+"_dadi.txt",dadi_vals_list=all_vals_recent,out_dir=param["out_dir"], title =p, xlim = False, ylim =False,name_pop=p)

#GQ distribution
if args.gq_distrib:
    if not os.path.exists(param["out_dir_gq_distrib"]):
        os.makedirs(param["out_dir_gq_distrib"])
    GQ_dict = res_pars[1]
    for p in param["name_pop"]:
        plot_distrib_gq(popid = p, gq = GQ_dict[p], out_dir_gq = param["out_dir_gq_distrib"] )

##SMCPP
if args.smcpp:
    if not os.path.exists(param["out_dir_smcpp"]):
        os.makedirs(param["out_dir_smcpp"])
    for p in param["name_pop"]:
        contigs = res_pars[2]
        smcpp(contigs = contigs, popid = p, pop_ind = param[p], vcf = param["vcf"], \
               out_dir = param["out_dir_smcpp"], mu = param["mut_rate"], gen_time = param["gen_time"], length_cutoff = param["length_cutoff"])
if args.plot_smcpp:
    for p in param["name_pop"]:
        plot_smcpp(popid = p, summary_file = param["plot_file_smcpp_"+p], out_dir = param["final_out_dir"])


##Gplot
if args.Gplot:
    for p in param["name_pop"]:
        Gplot(T_scaled_gen=param["out_dir_dadi"]+"popt_"+p+"_dadi.txt",gen_time=param["gen_time"],dadi_vals_list=dadi_output_parse(param["out_dir_dadi"]+"output_"+p+".dadi"),out_dir=param["out_dir"], title =p,name_pop=p,popid = p, summary_file2 = param["plot_file_smcpp_"+p],summary_file = param["summary_file_stw_"+p])



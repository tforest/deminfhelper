#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot  the plots
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import pandas as pd

def plot_sfs(sfs, plot_title, output_file):
    plt.bar(np.arange(1, len(sfs)+1), sfs)
    plt.title(plot_title)
    plt.xlabel('# of alternate alleles')
    plt.ylabel('# sites')
    plt.xticks(np.arange(1, len(sfs)+1, 1.0))
    plt.savefig(output_file)
    plt.close()


def plot_stairwayplot2(popid, summary_file, out_dir):
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
    plt.savefig(out_dir+popid+"_stw_inference.png")
    plt.close()

def plot_msmc2(popid, summary_file, mu, gen_time, out_dir):
    #scaled_left_bound = np.array()
    mu = float(mu)
    gen_time = float(gen_time)
    df = pd.read_csv(summary_file, delimiter='\t')
    # Scale Time
    # from units of the per-generation mutation rate to generations
    # to do this, divid by mutation rate
    scaled_left_bound = df['left_time_boundary']/mu * gen_time
    scaled_right_bound = df['right_time_boundary']/mu * gen_time
    # Pop Size
    scaled_pop = 1/df['lambda']
    pop_size_change = scaled_pop / (2*mu)
    plt.plot(scaled_left_bound, pop_size_change,color="blue")
    plt.title(popid)
    plt.xlabel('Time (in years)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+"_msmc2_inference.png")
    plt.close()

def plot_distrib_gq(popid, gq, out_dir_gq):
    gq = OrderedDict(sorted(gq.items()))
    names = list(gq.keys())
    values = list(gq.values())
    plt.figure(figsize = [10, 5])
    plt.bar(range(len(gq)), values, tick_label=names)
    plt.title(popid + " GQ distribution")
    plt.xlabel('GQ')
    plt.ylabel('Numbre of sites')
    plt.savefig(out_dir_gq+popid + "_gq_distrib.png")
    plt.close()


def plot_smcpp(popid, summary_file, out_dir):
    Ne=[]
    T=[]
    with open(summary_file) as input_file:
        line = input_file.readline()
        line = input_file.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T:
            Ne.append(float(line.split(',')[2]))
            T.append(float(line.split(',')[1]))
            line = input_file.readline()
    plt.plot(T,Ne,color="blue")
    plt.title(popid)
    plt.xlabel('Time (in years)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+"_smcpp_inference.png")
    plt.close()


def print_dadi_output_two_epochs(T_scaled_gen,dadi_vals_list,name_pop,out_dir, xlim, ylim, max_v = -10**6, title="Dadi pop. estimates"):
    with open(T_scaled_gen) as fo:
        line = fo.readline()
        fo.close()
    t_scaled_gen= line[1:-1]
    T_scaled_years=float(str.split(t_scaled_gen,", ")[-1])
    print("T_scaled_years",T_scaled_years)
    for dadi_params in dadi_vals_list:
        dadi_params[2][2]*=T_scaled_years
        dadi_params[2][3]*=T_scaled_years
    best_model = None
    #bad models :
    for elem in dadi_vals_list:
        popt = elem[2]
        if elem[1] > max_v:
            max_v = elem[1]
            best_model = elem[2]
        x =  [0, popt[3], popt[3], popt[3]+popt[2], popt[3]+popt[2], (popt[3]+popt[2])+0.01]
        y = [popt[1], popt[1], popt[0], popt[0], 1, 1]
        plt.plot(x, y, '--', alpha = 0.4)
    #best model :
    best_x = [0, best_model[3], best_model[3], best_model[3]+best_model[2], best_model[3]+best_model[2], (best_model[3]+best_model[2])+0.01]
    print("BEST",best_x)
    best_y = [best_model[1], best_model[1], best_model[0], best_model[0], 1, 1]
    plt.plot(best_x, best_y, 'r-', lw=2)
    plt.ylabel("Individuals (Na)")
    plt.xlabel("Time (years)")
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.title(title)
    plt.savefig(out_dir+"/inferences/output_plot_"+name_pop+"dadi.png")
    plt.close()

def Gplot(T_scaled_gen,gen_time,dadi_vals_list,name_pop,out_dir,popid, summary_file,summary_file2, max_v = -10**6, title="estimates",):
    import dadi
    with open(T_scaled_gen) as fo:
        line = fo.readline()
        fo.close()
    t_scaled_gen= line[1:-1]
    T_scaled_years=float(str.split(t_scaled_gen,", ")[-1])

    for dadi_params in dadi_vals_list:
        dadi_params[2][2]*=T_scaled_years
        dadi_params[2][3]*=T_scaled_years

    best_model = None
    for elem in dadi_vals_list:
        popt = elem[2]
        if elem[1] > max_v:
            max_v = elem[1]
            best_model = elem[2]
    best_x = [0, best_model[3], best_model[3], best_model[3]+best_model[2], best_model[3]+best_model[2], (best_model[3]+best_model[2])+0.01]
    best_y = [best_model[1], best_model[1], best_model[0], best_model[0], 1, 1]
    #with open(T_scaled_gen) as fo:
    #    line = fo.readline()
    #    fo.close()
    #t_scaled_gen= line[1:-1]
    #t_scaled_gen=str.split(t_scaled_gen,", ")[-1]
    Nanc=str.split(t_scaled_gen,", ")[0]
    #t_scaled_gen=[i*float(t_scaled_gen) for i in best_x]
    plt.plot(best_x, best_y,color="green")
    Ne=[]
    T=[]
    with open(summary_file2) as input_file2:
        line = input_file2.readline()
        line = input_file2.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T:
            Ne.append(float(line.split(',')[2]))
            T.append(float(line.split(',')[1]))
            line = input_file2.readline()
    plt.plot(T,[i/(2*Ne[-1]) for i in Ne],color="blue")

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
    plt.plot(T,[i/(2*Ne_med[-1]) for i in Ne_med],color="red")
    #plt.plot(T,[i/(2*Ne1[-1]) for i in Ne1],color="grey")
    #plt.plot(T,[i/(2*Ne3[-1]) for i in Ne3],color="grey")
    plt.title(popid)
    #plt.xlim(0,18000)
    plt.xlabel('Time (in years)')
    plt.ylabel('Na')
    plt.savefig(out_dir+popid+"_GPLOT_inference.png")
    plt.close()

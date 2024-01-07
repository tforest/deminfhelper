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
    
def plot_dadi_output_three_epochs(dadi_vals_list,name_pop,out_dir, mu, L, gen_time,
                                  xlim = None, ylim = None,
                                  max_v = -10**6, nb_plots_max = 10, title="Dadi pop. estimates"):
    best_model = None
    scenarios = {}
    for elem in dadi_vals_list:
        # elem is structured like this:
        # [iter_number, log_likelihood, [nuB, nuF, TB, TF], theta]
        log_lik = elem[1]
        # store scenarios sorted by their log_lik
        scenarios[float(log_lik)] = elem
    best_scenarios = sorted(scenarios.keys(), reverse = True)[:nb_plots_max]
    lines_x = []
    lines_y = []
    for scenario in best_scenarios:
        elem = scenarios[scenario]
        log_lik = elem[1]
        #popt = np.array(elem[2])
        theta = elem[3]
        Nanc = theta / (4*mu*L)
        # rescale to Nancestral effective size
        nuB = np.array(elem[2][0]) * Nanc
        nuF = np.array(elem[2][1]) * Nanc
        # Convert to years
        TB = np.array(elem[2][2]) * gen_time * Nanc
        TF = np.array(elem[2][3]) * gen_time * Nanc
        # now create x and y
        lines_x.append([0, TF, TF, TF+TB, TF+TB, (TF+TB)+0.01])
        lines_y.append([nuF, nuF, nuB, nuB, 1, 1])
    for i in range(1, len(lines_x)):
        x = lines_x[i]
        y = lines_y[i]
        plt.plot(x, y, '--', alpha = 0.4)
    #best model :
    best_x = lines_x[0]
    best_y = lines_y[0]
    plt.plot(best_x, best_y, 'r-', lw=2)
    plt.ylabel("Individuals (Na)")
    plt.xlabel("Time (years)")
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.title(title)
    plt.savefig(out_dir+"output_plot_"+name_pop+"dadi.png")
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

def plot_pca(plink_eigenvec, plink_eigenval):
    # Load eigenvectors and eigenvalues
    with open(plink_eigenvec) as input:
        for line in input:
            if line.startswith("#"):
                continue
            else:
                num_cols = len(line.split())
                break
    eigenvectors = np.genfromtxt(plink_eigenvec, skip_header=1, usecols=range(2, num_cols))  # Skip header and first two columns
    eigenvalues = np.loadtxt(plink_eigenval)

    # Determine the number of components
    num_components = eigenvectors.shape[1]

    # Calculate the percentage of variance explained by each PC
    variance_explained = (eigenvalues / np.sum(eigenvalues)) * 100

    # Plot the first two principal components
    plt.figure(figsize=(8, 6))
    plt.scatter(eigenvectors[:, 0], eigenvectors[:, 1], c='blue', marker='o', s=50)
    plt.title(f'PCA: PC1 vs PC2 ({num_components} components)')
    plt.xlabel(f'PC1 ({variance_explained[0]:.2f}%)')
    plt.ylabel(f'PC2 ({variance_explained[1]:.2f}%)')
    plt.show()

    # Bar plot of explained variance
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, num_components + 1), variance_explained * 100, color='blue', alpha=0.7)
    plt.xlabel('Number of Components (K)')
    plt.ylabel('Proportion of Explained Variance (%)')
    plt.title(f'Explained Variance by Components ({num_components} components)')
    plt.show()

# -*- coding: utf-8 -*-
"""
parsing() : parses the vcd, compute the sfs if SFS=True, output the msmc input file
if MSMC=False
"""

import gzip
# from inferences import *

if __package__ is None or __package__ == '':
    import inferences
    import sfs
else:
    from . import inferences
    from . import sfs

import re
from tqdm import tqdm  # Import tqdm for the progress bar

def parse_sfs(sfs_file):
    try:
        with open(sfs_file, 'r') as file:
            # Read the first line which contains information about the file
            num_individuals, mode, species_name = file.readline().strip().split()
            num_individuals = int(num_individuals)

            # print(f"Number of individuals: {num_individuals}")
            # print(f"Mode: {mode}")
            # print(f"Species name: {species_name}")

            # Read the spectrum data
            spectrum_data = list(map(int, file.readline().strip().split()))

            # Check if the number of bins in the spectrum matches the expected number
            if len(spectrum_data) != num_individuals:
                raise ValueError("Error: Number of bins in the spectrum doesn't match the expected number of individuals.")

            # print(f"Spectrum data: {spectrum_data}")

            # Read the mask data
            mask_data = list(map(int, file.readline().strip().split()))

            # Check if the size of the mask matches the number of bins in the spectrum
            if len(mask_data) != num_individuals:
                raise ValueError("Error: Size of the mask doesn't match the number of bins in the spectrum.")

            # print(f"Mask data: {mask_data}") 

            # Apply the mask to the spectrum
            masked_spectrum = [spectrum_data[i] for i in range(num_individuals) if not mask_data[i]]
            # print(f"Masked spectrum: {masked_spectrum}")

    except FileNotFoundError:
        print(f"Error: File not found - {file_path}")
    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"Error: {e}")
    # final return of SFS as a list
    return masked_spectrum
def get_contigs_lengths(param, length_cutoff=100000):
    contigs = []
    if 'length_cutoff' not in param.keys():
        param["length_cutoff"] = length_cutoff
    else:
        length_cutoff = param["length_cutoff"]

    with gzip.open(param["vcf"], 'rt') as vcf:
        line = vcf.readline()
        #print(line)
        while line != "":
            if line[0:8] == "##contig":
                # keep only contigs that are longer than the length_cutoff parameter
                contig_length = int(re.split('[=,]', line)[-1][:-2])
                if contig_length >= length_cutoff:
                    contigs.append(re.split('[=,]', line)[2])
            if line.startswith("#"):
                pass
            else:
                break
            line = vcf.readline()
    return contigs

def parsing(PARAM, SFS = False, GQ = False, SMCPP = False):
    # cutoff is the minimum size of each contig to be used
    # required for SMC++, as it works for contigs > 100kb or 1mb
    length_cutoff = int(PARAM["length_cutoff"])
    with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
        SFS_dict = {}
        GQ_dict = {}
        contigs = []
        Nseq = 0
        All_snp_count = 0
        Kept_snp_count = 0
        if SFS:
            # we initialize a sfs for each population
            for p in PARAM["name_pop"]:
                SFS_dict[p] = sfs.build_sfs(n=PARAM["n_"+str(p)], folded=PARAM["folded"], sfs_ini=True)
        if GQ:
            for p in PARAM["name_pop"]:
                        GQ_dict[p] = {}
        if SMCPP:
            contigs = []
        line = vcf.readline()
        pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar
        # we read all the lines of the vcf
        while line != "":
            if line[0:8] == "##contig":
                # keep only contigs that are longer than the length_cutoff parameter
                contig_length = int(re.split('[=,]', line)[-1][:-2])
                if contig_length >= length_cutoff:
                    contigs.append(re.split('[=,]', line)[2])
                    Nseq += contig_length
            if line[0:6] == "#CHROM":
                for p in PARAM["name_pop"]:
                    pos_ind_pop = []
                    for ind in PARAM[p]:
                        pos_ind_pop.append(line[:-1].split("\t").index(ind))
                    PARAM["pos_"+p] = pos_ind_pop
            if line.startswith("#"):
                pass
            if SFS or GQ:
                All_snp_count += 1
                if line[0] != "#" and ".:" not in line and "/." not in line and "./" not in line and ".|" not in line and "|." not in line and "," not in line.split("\t")[4]:    #we only keep the bi-allelique sites
                    Kept_snp_count += 1
                    split_line = line.split("\t")
                    if SFS:
                        for p in PARAM["name_pop"]:
                            SFS_dict[p] = sfs.build_sfs(n=PARAM["n_"+p], folded=PARAM["folded"],  sfs_ini = False, \
                                    line = split_line, sfs = SFS_dict[p], pos_ind = PARAM["pos_"+p])
                    if GQ:
                        for p in PARAM["name_pop"]:
                            GQ_dict[p] = distrib_GQ(GQ_pop = GQ_dict[p], line = split_line, pos_ind = PARAM["pos_"+p])
            line = vcf.readline()
            pbar.update(1)
    pbar.close()  # Close the progress bar when done
    L = (All_snp_count - Kept_snp_count) / All_snp_count * Nseq
    return SFS_dict, GQ_dict, contigs, round(L)

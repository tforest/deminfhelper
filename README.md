# DemInfHelper

## Overview

DemInfHelper is a demographic inference toolkit, specifically focused on population size change estimation.

The primary goal of this command line utility is to facilitate the estimation of population size changes from genomic variants. Users can choose from a suite of popular inference tools. Its modular design allows customization and parameter adjustments based on specific analysis needs and data specificity.

## Integrated Tools

1. **PSMC (Pairwise Sequentially Markovian Coalescent)**
   - **Description:** PSMC infers population size history from diploid sequences using a pairwise sequentially Markovian coalescent model. It explores scaled mutation and recombination rates, providing insights into demographic changes.
   - **Source:** [GitHub - lh3/psmc](https://github.com/lh3/psmc)
   - **Citation:** Li, H., and R. Durbin. 2011. "Inference of human population history from individual whole-genome sequences." Nature 475: 493–496. 

2. **MSMC2 (Multiple Sequentially Markovian Coalescent 2)**
   - **Description:** MSMC2 extends the MSMC model for inferring population size history and separation from whole-genome sequencing data. It provides accurate estimations for a large number of haplotypes.
   - **Source:** [GitHub - stschiff/msmc2](https://github.com/stschiff/msmc2)
   - **Citation:** Schiffels, S., and K. Wang. 2020. "MSMC and MSMC2: The Multiple Sequentially Markovian Coalescent." In Statistical Population Genomics, edited by J. Y. Dutheil, 147–166. Methods in Molecular Biology, Springer US, New York, NY.

3. **Stairway Plot v2**
   - **Description:** Stairway Plot v2 infers detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It can use both folded and unfolded SFSs and controls for overfitting.
   - **Source:** [GitHub - xiaoming-liu/stairway-plot-v2](https://github.com/xiaoming-liu/stairway-plot-v2)
   - **Citation:** Liu, X., and Y.-X. Fu, 2020 Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biology 21: 280.

4. **dadi-cli (Flexible Python Package for Demographic Inference)**
   - **Description:** dadi-cli provides a user-friendly command-line interface for dadi, a flexible Python package for inferring demographic history and the distribution of fitness effects from population genomic data.
   - **Source:** [dadi-cli documentation](https://dadi-cli.readthedocs.io)
   - **Citation:** Gutenkunst, R. N., R. D. Hernandez, S. H. Williamson, and C. D. Bustamante. 2009. "Inferring the Joint Demographic History of Multiple Populations from Multidimensional SNP Frequency Data." PLOS Genetics 5: e1000695.
 
5. **SMC++ (Sequentially Markovian Coalescent++)**
   - **Description:** SMC++ is a tool for estimating the size history of populations from whole-genome sequence data. It offers various subcommands for data conversion, estimation, cross-validation, and joint demography analysis.
   - **Source:** [GitHub - popgenmethods/smcpp](https://github.com/popgenmethods/smcpp)
   - **Citation:** Terhorst, J., J. A. Kamm, and Y. S. Song. 2017. "Robust and scalable inference of population history from hundreds of unphased whole genomes." Nat Genet 49: 303–309.

## Repository Structure

- **bin/:** Contains executable files for the demographic inference tools, essentially the StairwayPlot2 executable, as it is not available in Conda.
- **scripts/:** Contains external`generate_multihetsep.py` from the `msmc-tools` repository.
- **deminfhelper.py:** Main script. 
- **environment.yml:** Conda environment file specifying dependencies to run the tools.
- **inferences.py**: Definitions of wrapers to run all the inference tools.
- **parsing.py**: Parsing functions for config or input files.
- **plots.py**: Plotting functions.
- **sfs.py**: Functions to produce the Site Frequency Spectrum (SFS) from VCF files.

## Usage

1. **Installation:**
   - Get DemInfHelper. You can clone the repository straight from Github to get the latest version :
     ```bash
     git clone https://github.com/tforest/deminfhelper.git
     ```

2. **Dependencies:**

It is recommended to get all the required dependencies using Conda, inside a dedicated conda environment, using the provided environment.yml file:
   ```bash
     conda env create -f environment.yml
     conda activate deminfhelper
  ``` 


3. **Configuration:**
   - Create a configuration in a yaml file.
 
 Example **cfg.yml** file:
 ```yaml
    ### CONFIGURATION FILE for DemInfHelper
    ## Output top level directory
    out_dir: ./picus_viridis/
    ## Input data
    vcf: ./VCF/picus_viridis_raw.vcf.gz
    ref_genome: ./GENOME/GCA_033816785.1.fna
    ## L: effective sequence length genotyped
    L: 1244800746
    ## Population(s)
    name_pop: picus_viridis
    npop: 1
    # POP 1
    # Note: If you don't specify samples, it will use all samples in the VCF
    picus_viridis: SAMN38508692, SAMN38508693, SAMN38508694, SAMN38508695, SAMN38508696, SAMN38508697, SAMN38508698, SAMN38508699, SAMN38508701, SAMN38508702
    # PARAM
    gen_time: 2
    mut_rate: 3e-8
    ## SFS
    out_dir_sfs: ./picus_viridis/
    folded: True
    # Only keep contigs >=10kb for the analyses
    length_cutoff: 10000
    ## SFS if previously computed
    path_to_sfs: ./picus_viridis/SFS_picus_viridis.fs
    ## StairwayPlot2
    path_to_stairwayplot2: ../deminfhelper/bin/stairway_plot_es/
    blueprint_template: ../deminfhelper/bin/template.blueprint
    out_dir_stairwayplot2: ./picus_viridis/output_stairwayplot2/
    summary_file_stw: ./picus_viridis/output_stairwayplot2/picus_viridis/picus_viridis.final.summary
    ## SMC++
    out_dir_smcpp: ./picus_viridis/output_smcpp/
    plot_file_smcpp: ./picus_viridis/output_smcpp/picus_viridis_inference.csv
    # MSMC2
    msmc2_kwargs: -i 25 -p 1*2+25*1+1*2+1*3
    out_dir_msmc2: ./picus_viridis/output_msmc2/
    ## PSMC
    out_dir_psmc: ./picus_viridis/output_psmc/
    psmc_kwargs: -N25 -t15 -r5 -p "4+25*2+4+6"
    plot_psmc_kwargs: -x 10**4 -X 150000    
    ## dadi parameters
    out_dir_dadi: ./picus_viridis/output_dadi/
    optimizations: 100
    lower_bound: 1, 1, 0.05, 0.01
    p0: 0.01, 0.001, 0.01, 0.01
    upper_bound: 10, 4, 0.1, 10
    ## GQ distribution
    out_dir_gq_distrib: ./picus_viridis/output_stats/
    ## FINAL INFERENCES
    final_out_dir: ./picus_viridis/inferences/
    ## Stats
    out_dir_stats: ./picus_viridis/output_stats/
    ## PCA K-means clustering
    n_clust_kmeans: 3
    ## Filtering contig based on their name using regex.
    # Filter contigs starting by "CM":
    contig_filter: CM.*
    # Note: To omit filtering, use .*
 ```
4. **Usage:**
```bash
$ ./deminfhelper.py -h

usage: deminfhelper.py [-h] [--config_file CONFIG_FILE] [--cpus CPUS] [--sfs] [--sfs_transformed] 
[--plot_sfs] [--stairwayplot2] [--plot_stairwayplot2] [--dadi]
                       [--plot_dadi] [--msmc2] [--plot_msmc2] 
                       [--msmc2_kwargs MSMC2_KWARGS] [--psmc] [--plot_psmc] 
                       [--plot_psmc_kwargs PLOT_PSMC_KWARGS]
                       [--psmc_kwargs PSMC_KWARGS] [--gq_distrib] [--smcpp] 
                       [--plot_smcpp] [--Gplot] [--folded] [--pca] [--plot_pca] 
                       [--n_clust_kmeans N_CLUST_KMEANS]
                       [--popid POPID] [--samples SAMPLES] [--vcf VCF] [--gentime GENTIME] 
                       [--mu MU] [--out OUT] [--out_dir_stats OUT_DIR_STATS] [--L L] [--n N]
                       [--contig_filter CONTIG_FILTER]

Computes the sfs from the vcf and runs demography inference softwares.

optional arguments:
  -h, --help            show this help message and exit
  --config_file CONFIG_FILE
                        path to the configuration file
  --cpus CPUS           number of CPU threads to use
  --mask MASK           Keep only regions specified in a given BED file.
  --sfs                 to compute the sfs
  --sfs_transformed     to normalize the sfs
  --plot_sfs            to plot the sfs
  --stairwayplot2       to run stairwayplot2
  --plot_stairwayplot2  to run stairwayplot2
  --dadi                to run dadi: the sfs must not be transformed
  --plot_dadi           to create popsize plot from dadi output.
  --msmc2               to run msmc2
  --plot_msmc2          to plot msmc2
  --msmc2_kwargs MSMC2_KWARGS
                        Optional args for MSMC2 (list separated by ';', 
                        eg. --msmc2_kwarg -i 25; -p 1*2+25*1+1*2+1*3 and so on).
  --psmc                to run PSMC
  --plot_psmc           to plot psmc inference
  --plot_psmc_kwargs PLOT_PSMC_KWARGS
                        Params to pass to plot_psmc.pl script. Precise them inside quotes 
                        (eg. '-x 10**4 -X 150000')
  --psmc_kwargs PSMC_KWARGS
                        PSMC params (eg. --psmc_kwarg -p 10+22*2+4+6; -N30 -t15 -r5 and so on)
  --gq_distrib          to compute the GQ (genotype quality) distribution
  --smcpp               run smcpp
  --plot_smcpp          to plot smcpp inference
  --Gplot               to plot all inferences on the same graph
  --folded              Fold the SFS. Default: True
  --pca                 Compute PCA using Plink2
  --plot_pca            Compute PCA using Plink2
  --n_clust_kmeans N_CLUST_KMEANS
                        Defines the number of k clusters for the k-means on the PCA
  --popid POPID         a population identifier; eg. species name
  --samples SAMPLES     a list of samples to use in the VCF. By default all samples are taken.
  --vcf VCF             the vcf file path, only gz format
  --gentime GENTIME     the generation time of the species in years. Eg 1.5
  --mu MU               the average mutation rate per site per generation
  --out OUT             Output path of the analysis. Will create a dir.
  --out_dir_stats OUT_DIR_STATS
                        Output path of the stats. Will create a dir.
  --L L                 The actual size of the genotyped sequence length that produced the VCF, 
                        before any filters.
  --n N                 Number of sequences that produced the VCF. Eg for 5 dipl individuals, n=10.
  --contig_filter CONTIG_FILTER
                        Keep only contigs satisfiying that regular expression.

```

5. **Example use:**
- Compute the SFS
     ```bash
     deminfhelper.py --config_file cfg.yml --sfs
     ```
- Plot the SFS
     ```bash
     deminfhelper.py --config_file cfg.yml --plot_sfs
     ```
- Plot the Genotype Quality (GQ) distribution
     ```bash
     deminfhelper.py --config_file cfg.yml --gq_distrib
     ```
- Run inference(s)
     ```bash
     deminfhelper.py --config_file cfg.yml --stairwayplot2
     deminfhelper.py --config_file cfg.yml --dadi
     deminfhelper.py --config_file cfg.yml --smcpp
     deminfhelper.py --config_file cfg.yml --msmc2 --cpus 8
     deminfhelper.py --config_file cfg.yml --psmc --cpus 8
     ```
 - Plot the inference(s)
      ```bash
     deminfhelper.py --config_file cfg.yml --plot_stairwayplot2
     deminfhelper.py --config_file cfg.yml --plot_dadi
     deminfhelper.py --config_file cfg.yml --plot_smcpp
     deminfhelper.py --config_file cfg.yml --plot_msmc2 --cpus 8
     deminfhelper.py --config_file cfg.yml --plot_psmc --cpus 8
     ```


## Contribute

If you encounter issues or have suggestions, please open an issue on the [GitHub repository](https://github.com/tforest/deminfhelper/issues).


## Citation

If you use this tool, please consider citing the following publications:

- Gutenkunst, R. N., R. D. Hernandez, S. H. Williamson, and C. D. Bustamante. 2009. "Inferring the Joint Demographic History of Multiple Populations from Multidimensional SNP Frequency Data." PLOS Genetics 5: e1000695.
- Li, H., and R. Durbin. 2011. "Inference of human population history from individual whole-genome sequences." Nature 475: 493–496.
- Schiffels, S., and K. Wang. 2020. "MSMC and MSMC2: The Multiple Sequentially Markovian Coalescent." In Statistical Population Genomics, edited by J. Y. Dutheil, 147–166. Methods in Molecular Biology, Springer US, New York, NY.
- Terhorst, J., J. A. Kamm, and Y. S. Song. 2017. "Robust and scalable inference of population history from hundreds of unphased whole genomes." Nat Genet 49: 303–309.
- Liu, X., and Y.-X. Fu, 2020 Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biology 21: 280.

# Code and Script for CAPS paper

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4321048.svg)](https://doi.org/10.5281/zenodo.4321048)

Authors: Yibin Liu,<sup>1,2,5*</sup> Zhiyuan Hu,<sup>3,4,*</sup> Jingfei Cheng,<sup>1,2</sup> Paulina Siejka-Zielińska,<sup>1,2</sup> Jinfeng Chen,<sup>1,2</sup> Masato Inoue,<sup>1,2</sup> Ahmed Ashour Ahmed,<sup>3,4</sup> Chun-Xiao Song<sup>1,2,†</sup>  
Affiliations:  
1 Ludwig Institute for Cancer Research, Nufﬁeld Department of Medicine, University of Oxford, Oxford OX3 7FZ, UK  
2 Target Discovery Institute, Nufﬁeld Department of Medicine, University of Oxford, Oxford OX3 7FZ, UK  
3 Ovarian Cancer Cell Laboratory, MRC Weatherall Institute of Molecular Medicine, University of Oxford, Oxford OX3 9DS, UK  
4 Nuffield Department of Women’s & Reproductive Health, University of Oxford, Oxford OX3 9DU, UK  
5 Present address: Exact Sciences Innovation, Innovation Building, Oxford OX3 7FZ, UK  
*These authors contributed equally to this work.  
†Corresponding author. Email: chunxiao.song@ludwig.ox.ac.uk  


## Citation

## Folder Structure

**snake_pipeline/**
* Snakemake pipeline to preprocessing the fastq files and generate methylation calling files
* The development version asTair-3.3.1 of astair was installed by `python -m pip install --user astair==3.3.1`.

**analysis/**
* Statistical analysis and visualisation

## Data Preprocessing: Snakemake Workflow

**Description:**
* Preprocessing TAPS&beta;, CAPS, PS and PS-c. 
* Reading in fastq files and generated methylation calling results.
* `asTair vesion: 3.3.1` used

**Files:**
* snake_pipeline/

**Steps included:**
* Trim reads
* Alignment (bwa mem)
* Call methylation sites on spikein
* Deduplicate
* Clip overlap
* Call methylation sites on mm9 genome (SNV excl.)
* Mask artifact-prone regions
* [QC] Phred visualisation
* [QC] Counting reads before and after trimming, alignment and deduplication


## Downstream Analysis

The following analysis was aimed to summarise and visualise the methylation calling results. I listed the scripts correpsonding to the analysis.

### Spike-in analysis

1. analysis/spikein_analysis.Rmd (local)

### Visualise chemical results (related to Figure S3a, b)

1. analysis/caps_HPLC_boxplot.R: plot the HPLC results


### TAPS&beta; analysis (related to Figure 1d, e)

1. analysis/tapsbeta_5mc_benchmarking_taps_oxbs.R (remote)
1. analysis/tapsbeta_5mc_benchmark_runr.sh (remote): to run the above Rscript


### CAPS: benchmarking analysis (related to Figure 2e,f)

1. analysis/mlml_tapsbeta_preprocessing.R
1. analysis/mlml.sh (run preprocessing R script and run MLML)
1. analysis/caps_analysis_preprocessing.R (ran by caps_analysis_bin_chrs.sh): preprocessing CAPS data, mlml results, ACE data and TABseq data.
1. analysis/caps_analysis_bin_chrs.sh: bedtools map to count Ts and Cs in 10-kb bins
1. analysis/caps_analysis_benchmarking.R (local)

### CAPS: CpG islands coverage analysis (related to Figure S6)

1. analysis/CpG_island_prepare_bins.R
1. analysis/CpG_island_analysis.sh
1. analysis/CpG_island_analysis.R


### Subtraction analysis (related to Figure 3a)

1. analysis/subtraction_preprocessing.R 
1. analysis/subtraction_mlml.sh
1. analysis/subtraction_diagram.R (local)


### IGV visualisation (related to Figure 3b, S7)

1. analysis/igv_mod2bed.R: astair output to bedgraph
1. analysis/igv_bed2bw.sh: bedgraph to bigwig


### CAPS: Colocolization analysis with genomic regulatory elements (related to Figure 3c, d)

1. analysis/caps_genomic_element_analysis.R - part 1 (commanded out)
1. analysis/caps_genomic_element_analysis.sh (remotely)
1. analysis/caps_genomic_element_analysis.R - part 2

### PS analysis (related to Figure 4f,4g,S8,9)

1. ps_broadpeak_final.sh: methylation signal around histone modification peaks
1. ps_hmm_dts_final.sh: methylation signal in predicted genomic elements
1. ps_astair2bw.sh: for IGV visualisation


## Session info

### Server enviroment

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.10 (Final)

Matrix products: default
BLAS:   /usr/lib64/libblas.so.3.2.1
LAPACK: /usr/lib64/atlas/liblapack.so.3.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] ggplot2_3.2.1     viridis_0.5.1     viridisLite_0.3.0 MASS_7.3-51.5    
[5] data.table_1.12.6

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.3       withr_2.1.2      assertthat_0.2.1 dplyr_0.8.4     
 [5] crayon_1.3.4     grid_3.6.0       R6_2.4.1         lifecycle_0.1.0 
 [9] gtable_0.3.0     magrittr_1.5     scales_1.1.0     pillar_1.4.3    
[13] rlang_0.4.1      lazyeval_0.2.2   glue_1.3.1       purrr_0.3.3     
[17] munsell_0.5.0    compiler_3.6.0   pkgconfig_2.0.3  colorspace_1.4-1
[21] tidyselect_0.2.5 gridExtra_2.3    tibble_2.1.3    
```

### Local enviroment

```
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] data.table_1.12.8

loaded via a namespace (and not attached):
[1] compiler_3.6.2     tools_3.6.2        yaml_2.2.1         KernSmooth_2.23-16
```


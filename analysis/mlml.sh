#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq


## Zhiyuan Hu
## May 2020
## updated 26 july 2020

## mlml of TAPS and TAPSbeta
module load R/3.6.0-newgcc
Rscript mlml_tapsbeta_preprocessing.R

cd $HOME/projects/taps/analysis_caps/tidy_data/mlml/
gzip -d tidy_tapsbeta_mlml.tsv.gz
gzip -d tidy_taps_mlml.tsv.gz

mlml -v -u tidy_taps_mlml.tsv -m tidy_tapsbeta_mlml.tsv -o taps_tapsbeta_mlml_output.tsv

gzip taps_tapsbeta_mlml_output.tsv
gzip tidy_tapsbeta_mlml.tsv
gzip tidy_taps_mlml.tsv

## preprocessing BS-seq file
cd $HOME/projects/taps/benchmark/oxbs
zcat GSM4194641_Control_gDNA_100ng.cov.gz | sort -k 1,1 -k2,2n - | gzip > GSM4194641_Control_gDNA_100ng_sorted.cov.gz

## mlml of BS-seq and TAB-seq
mlml -v -u tidy_bs_chrall_mlml2.tsv -h tidy_tabseq_chrall_mlml2.tsv -o bs_tabseq_mlml_output_chrall.tsv
gzip bs_tabseq_mlml_output_chrall.tsv




## Using mlml from the MethPipe
## mlml parameters:
##  -o, -output     Name of output file (default: stdout) 
##  -u, -bsseq      Name of input BS-Seq methcounts file 
##  -h, -tabseq     Name of input TAB-Seq methcounts file 
##  -m, -oxbsseq    Name of input oxBS-Seq methcounts file 
##  -t, -tolerance  EM convergence threshold. Default 1e-10 
##  -a, -alpha      significance level of binomial test for each site. Default 
##                  0.05 
##  -v, -verbose    print run statistics 

## Input data format
## chr1    1       +       CpG     0.363636363636  11
## chr1    3       +       CpG     1.0     7

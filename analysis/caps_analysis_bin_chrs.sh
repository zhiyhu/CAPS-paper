#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq


## calculate 10-kb binned raw signals of 5hmC
## Zhiyuan Hu
## 24 Apr 2020
## last modified 26 july 2020

## output file: tidy_data/bins_out/5hmC_caps_mlml_ace_tab_mods_10kb_bin.bed.gz

module load bedtools/2.25.0
module load R/3.6.0-newgcc

## preprocessing data

echo 'Preprocessing data'
Rscript $HOME/projects/taps/analysis_caps/R/caps_analysis_preprocessing.R 2>&1 | tee > logs/caps_analysis_preprocessing_20200726.log

## count sum(C)/sum(C+T) in 10-kb bins

bin=$HOME/projects/taps/ref/genome.10kb.bed.gz
outpath=$HOME/projects/taps/analysis_caps/tidy_data/bins_out

in=$HOME/projects/taps/analysis_caps/tidy_data/5hmC_caps_mlml_ace_tab_mods.bed.gz
out=${outpath}/5hmC_caps_mlml_ace_tab_mods_10kb_bin.bed.gz
zcat $in | grep -v "strand" - | bedtools map -a $bin -b - -o sum -c 5,6,7,8,9,10,21,22,18,19,20 -null "NA" | gzip > $out


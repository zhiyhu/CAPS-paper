#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq

## annotate CpG with the genomic element annotation
## Zhiyuan Hu
## 20 May 2020

module load bedtools/2.25.0 

## the mESC_cStates_HMM.bed.gz was generated by R/genomic_element_analysis.R

## sort HMM file first (don't need to rerun this)
## bed=$HOME/projects/taps/ref/mESC_cStates_HMM.bed.gz
## out=$HOME/projects/taps/ref/mESC_cStates_HMM_sorted.bed
## zcat $bed | sort -k1,1 -k2,2n > $out
## gzip $out

## sort caps methy data
dir=$HOME/projects/taps/analysis_caps/methy_data
bed=${dir}/caps_CpG_mods.bed.gz
out=${dir}/caps_CpG_mods_sorted.bed
zcat $bed | sort -k1,1 -k2,2n > $out
gzip $out

## annoate caps methy data
dir=$HOME/projects/taps/analysis_caps/methy_data
bedb=$HOME/projects/taps/ref/mESC_cStates_HMM_sorted.bed.gz
beda=${dir}/caps_CpG_mods_sorted.bed.gz
out=${dir}/caps_CpG_mods_annotated.bed.gz
bedtools intersect -sorted -wo -a $beda -b $bedb | gzip > $out

## run R script

module load R/3.6.0-newgcc

Rscript $HOME/projects/taps/analysis_caps/R/caps_genomic_element_analysis.R 2>&1 | \
tee > $HOME/projects/taps/analysis_caps/R/logs/caps_genomic_element_analysis_20200726.log



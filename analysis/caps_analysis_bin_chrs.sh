#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq
## The following two lines will send an email notification when the job is
## Ended/Aborted/Suspended - Please replace "UserName" with your own username.
#$ -M zyhu
#$ -m aes

## calculate 10-kb binned raw signals of 5hmC
## Zhiyuan Hu
## 24 Apr 2020
## last modified 26 july 2020


module load bedtools/2.25.0
module load R/3.6.0-newgcc

## preprocessing data

echo 'Preprocessing data'
Rscript $HOME/projects/taps/analysis_caps/R/caps_analysis_preprocessing.R 2>&1 | tee > logs/caps_analysis_preprocessing_20200726.log

## count sum(C)/sum(C+T) in 10-kb bins

bin=/home/obgynae/zyhu/projects/taps/ref/genome.10kb.bed.gz
outpath=$HOME/projects/taps/analysis_caps/tidy_data/bins_out

in=$HOME/projects/taps/analysis_caps/tidy_data/5hmC_caps_mlml_ace_tab_mods.bed.gz
out=${outpath}/5hmC_caps_mlml_ace_tab_mods_10kb_bin.bed.gz
zcat $in | grep -v "strand" - | bedtools map -a $bin -b - -o sum -c 5,6,7,8,9,10,21,22,18,19,20 -null "NA" | gzip > $out









################################################################################################
################################################################################################
################################################################################################
### below is obsolete

in=$HOME/projects/taps/analysis_caps/tidy_data/5hmC_caps_ace_tab_mods.bed.gz
out=$outpath/5hmC_caps_ace_tab_mods_10kb_bin.bed
zcat $in | bedtools map -a $bin -b - -o sum -c 5,6,7,8,9,10,14,15,16 -null "NA" > $out
gzip $out

in=$HOME/projects/taps/analysis_caps/tidy_data/5hmC_tapsmlml_ace_tab_mods.bed.gz
out=$outpath/5hmC_tapsmlml_ace_tab_mods_10kb_bin.bed
zcat $in | bedtools map -a $bin -b - -o sum -c 5,6,7,8,9,10,13,14 -null "NA" > $out
gzip $out



## https://www.biostars.org/p/70577/

## chrom.sizes format
## chrI  15072421
## chrII 15279323
## ...
## chrX  17718854
## chrM  13794

## bedtools intersect https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

## bed format http://genome.ucsc.edu/FAQ/FAQformat#format1
## chrom	chromStart	chromEnd	name	score	strand



######################
## below are older codes
######################

module load bedtools/2.25.0

winsize=10000 # window size

## make bins
refpath=/home/obgynae/zyhu/projects/taps/ref
chromsizes=${refpath}/chrom_sizes.tsv
bin=${refpath}/genome.10kb.bed
# bedtools makewindows -g $chromsizes -w $winsize > $bin

## bedtools map
s=YL712508hm ## sample
dir=/home/obgynae/zyhu/projects/taps/caps/R/tidy_data
input=${dir}/${s}_dedup_mCtoT_chrAll_CpG_masked_columnSeleted.tsv ## tidy data

o1=${dir}/${s}_${winsize}bp_binned_sum_mod_reads.tsv
o2=${dir}/${s}_${winsize}bp_binned_sum_unmod_reads.tsv
o3=${dir}/${s}_${winsize}bp_binned_mean_mod_levels.tsv

bedtools map -a $bin -b $input -o sum -c 5 > $o1 # sum of mod reads
bedtools map -a $bin -b $input -o sum -c 6 > $o2 # sum of unmod reads
bedtools map -a $bin -b $input -o mean -c 4 > $o3 # mean of mod levels

## ACEseq
in=/home/obgynae/zyhu/projects/taps/caps/R/tidy_data/ACE/GSE116016_ACE-Seq_WT.ESC.mm9_CG.txt.gz
out=/home/obgynae/zyhu/projects/taps/caps/R/tidy_data/ACE_dedup_mCtoT_chrAll_CpG_masked_columnSeleted.tsv
zcat $in | awk 'BEGIN {FS="\t"; OFS="\t"}; {print $1,$2,$3,$4,$5 * $4,$5 - $5 * $4, $6}' > $out 


## CAPS
## ZHiyuan Hu
## 15 May 2020
module load bedtools/2.25.0

refpath=/home/obgynae/zyhu/projects/taps/ref
bin=${refpath}/genome.10kb.bed
in=$HOME/projects/taps/analysis_caps/tidy_data/5hmC_tapsmlml_ace_tab_mods.bed.gz
outpath=$HOME/projects/taps/analysis_caps/tidy_data/bins_out

##bl=$HOME/projects/taps/ref/mm9-blacklist.bed.gz # blacklist
##intersectBed -v -a $bl -b $bin > ${refpath}/genome.10kb_masked.bed


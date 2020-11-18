#!/bin/sh
# Astair output to bigwig by depth cutoff of 10
# Zhiyuan Hu
# 24 Oct 2020; last modified 18 Nov 2020

module load ucsctools/385
wd=$HOME/projects/taps/caps/astair_output
chrsize=$HOME/projects/taps/ref/chrom_sizes.tsv

# cutoff = 10
for smp in YL708508f YL709502c; do
infile=${wd}/mm9/${smp}_dedup_chrselected.clipOverlap_mCtoT_CpG_masked.mods.gz
bedfile=${wd}/bigwig/${smp}_CpG_modlevel_cutoff10.bedGraph
bwfile=${wd}/bigwig/${smp}_CpG_modlevel_cutoff10.bw

# filter by SNV and depth
# sort
# select conversion rate
# print to bedgraph
zcat $infile | \
awk -F"\t" '/No/ && ($5+$6)>=10' | \
cut -f 1,2,3,4 | \
sort -k1,1 -k2,2n > $bedfile

bedGraphToBigWig $bedfile $chrsize $bwfile
done
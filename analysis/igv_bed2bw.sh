#!/bin/sh

# bed to bigwig
# Output file: tidy_data/igv/*.bigwig

## refer to https://www.biostars.org/p/150036/

##########################
## TAPS, TAPSbeta and CAPS
## mod.gz to bedgraph (bed.gz)

# module load R/3.6.0-newgcc
# Rscript igv_mod2bed.R

####################
## grep chromosomes 4, 6 and 18
chr=chr4

cd $HOME/projects/taps/analysis_caps/tidy_data/igv

zcat taps_modlevel.bed.gz     | grep $chr | gzip > taps_modlevel_${chr}.bed.gz
zcat tapsbeta_modlevel.bed.gz | grep $chr | gzip > tapsbeta_modlevel_${chr}.bed.gz
zcat caps_modlevel.bed.gz     | grep $chr | gzip > caps_modlevel_${chr}.bed.gz
zcat ace_modlevel.bed.gz      | grep $chr | gzip > ace_modlevel_${chr}.bed.gz
zcat tabseq_modlevel.bed.gz   | grep $chr | gzip > tabseq_modlevel_${chr}.bed.gz

##########################
## bedgraph to wigwig
module load ucsctools/2.0

chrsize=/home/obgynae/zyhu/projects/taps/ref/chrom_sizes_sorted.tsv

for s in taps tapsbeta caps # ace tabseq
do
for tp in modlevel
do
fn=$HOME/projects/taps/analysis_caps/tidy_data/igv/${s}_${tp}_${chr}

gzip -d ${fn}.bed.gz 
bedGraphToBigWig ${fn}.bed $chrsize ${fn}.bigwig
gzip ${fn}.bed

done
done



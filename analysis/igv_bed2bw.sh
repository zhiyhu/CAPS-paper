#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq
## The following two lines will send an email notification when the job is
## Ended/Aborted/Suspended - Please replace "UserName" with your own username.
#$ -M zyhu
#$ -m aes

## refer to https://www.biostars.org/p/150036/

# awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' myBed.bed > myFile.bedgraph
# sort -k1,1 -k2,2n myFile.bedgraph > myFile_sorted.bedgraph

##########################
## TAPS, TAPSbeta and CAPS
## mod.gz to bedgraph (bed.gz)

# module load R/3.6.0-newgcc
# Rscript igv_mod2bed.R

####################
## grep chromosomes 4, 6 and 18
chr=chr4

cd /home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv

zcat taps_modlevel.bed.gz | grep $chr | gzip >  taps_modlevel_${chr}.bed.gz
# zcat taps_depth.bed.gz    | grep $chr | gzip >  taps_depth_${chr}.bed.gz

zcat tapsbeta_modlevel.bed.gz | grep $chr | gzip >  tapsbeta_modlevel_${chr}.bed.gz
# zcat tapsbeta_depth.bed.gz    | grep $chr | gzip >  tapsbeta_depth_${chr}.bed.gz

zcat caps_modlevel.bed.gz | grep $chr | gzip > caps_modlevel_${chr}.bed.gz
# zcat caps_depth.bed.gz    | grep $chr | gzip > caps_depth_${chr}.bed.gz

# zcat ace_modlevel.bed.gz | grep $chr | gzip > ace_modlevel_${chr}.bed.gz
# zcat ace_depth.bed.gz    | grep $chr | gzip > ace_depth_${chr}.bed.gz

# zcat tabseq_modlevel.bed.gz | grep $chr | gzip > tabseq_modlevel_${chr}.bed.gz
# zcat tabseq_depth.bed.gz    | grep $chr | gzip > tabseq_depth_${chr}.bed.gz


##########################
## bedgraph to wigwig
module load ucsctools/2.0

# chrsize=/home/obgynae/zyhu/projects/taps/ref/chrom_sizes.tsv
# sort -k1,1 $chrsize > /home/obgynae/zyhu/projects/taps/ref/chrom_sizes_sorted.tsv
chrsize=/home/obgynae/zyhu/projects/taps/ref/chrom_sizes_sorted.tsv

for s in taps tapsbeta caps # ace tabseq
do
for tp in modlevel
do
fn=/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/${s}_${tp}_${chr}

gzip -d ${fn}.bed.gz ## | sort -k1,1 -k2,2n  > ${fn}_sorted.bedgraph
bedGraphToBigWig ${fn}.bed $chrsize ${fn}.bigwig
gzip ${fn}.bed

done
done



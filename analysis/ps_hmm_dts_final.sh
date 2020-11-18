#!/bin/sh
# run ChromHMM annotations with deeptools
# Zhiyuan Hu
# 2 Nov 2020; last modified 2 Nov 2020

eval "$(/t1-data/data/ahmedgroup/zhiyuan/miniconda3/bin/conda shell.bash hook)"
conda init
conda activate dt

hmm=$HOME/projects/taps/ref/mESC_cStates_HMM_sorted.bed.gz
wd=$HOME/projects/taps/caps/deeptools
bl=$HOME/projects/taps/ref/mm9-blacklist.bed

smp3=YL708508f
smp4=YL709502c
depth=1
bs=100
bw1=$HOME/projects/taps/caps/astair_output/bigwig/${smp3}_CpG_modlevel_cutoff${depth}.bw
bw2=$HOME/projects/taps/caps/astair_output/bigwig/${smp4}_CpG_modlevel_cutoff${depth}.bw


for anno in Strong_Enhancer Active_Promoter Heterochrom Repressed; do

# not skip zeros
mat=${wd}/ps-psc/chromhmm/${smp3}.${smp4}_${anno}_cutoff${depth}_bin${bs}_blfiltered.tab.gz
pdf=${wd}/ps-psc/chromhmm/${smp3}.${smp4}_${anno}_cutoff${depth}_bin${bs}_blfiltered.pdf

zcat $hmm  | grep $anno > ${wd}/hmm/mESC_chromHMM_${anno}.bed

computeMatrix reference-point \
--referencePoint center -p 4 \
-a 3000 -b 3000 -bs $bs \
-bl $bl \
-S $bw1 $bw2 \
-R ${wd}/hmm/mESC_chromHMM_${anno}.bed \
-o $mat 

plotProfile \
-m $mat \
-out $pdf \
--samplesLabel $smp3 $smp4 \
--legendLocation upper-right \
--refPointLabel 0 \
--plotFileFormat pdf \
--perGroup \
--plotTitle $anno \
--yMin 0.0015 \
--yMax 0.0030

done

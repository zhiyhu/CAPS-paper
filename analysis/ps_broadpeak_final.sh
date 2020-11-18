#!/bin/sh
# deeptools analysis
# Zhiyuan Hu
# 24 Oct 2020; last modified 12 Nov 2020

eval "$(/t1-data/data/ahmedgroup/zhiyuan/miniconda3/bin/conda shell.bash hook)"
conda init
conda activate dt

smp3=YL708508f
smp4=YL709502c
wd=$HOME/projects/taps/caps/astair_output
bl=$HOME/projects/taps/ref/mm9-blacklist.bed
bedpath=$HOME/projects/taps/ref/encode
for hist in h3k4me1 h3k4me3; do
filetype=broadPeak 

bed=${bedpath}/e14_${hist}.${filetype}

depth=1
bs=100

bw1=${wd}/bigwig/${smp3}_CpG_modlevel_cutoff${depth}.bw
bw2=${wd}/bigwig/${smp4}_CpG_modlevel_cutoff${depth}.bw

# not skip zeros
mat=$HOME/projects/taps/caps/deeptools/ps-psc/final.${filetype}/${smp3}.${smp4}_${hist}.${filetype}_cutoff${depth}_bin${bs}.tab.gz
png=$HOME/projects/taps/caps/deeptools/ps-psc/final.${filetype}/${smp3}.${smp4}_${hist}.${filetype}_cutoff${depth}_bin${bs}.pdf

computeMatrix reference-point \
--referencePoint center -p 4 \
-a 3000 -b 3000 -bs $bs \
-bl $bl \
-S $bw1 $bw2 \
-R $bed \
-o $mat 

plotProfile \
-m $mat \
-out $png \
--samplesLabel $smp3 $smp4 \
--legendLocation upper-right \
--refPointLabel 0 \
--plotFileFormat pdf \
--perGroup \
--plotTitle $hist

done

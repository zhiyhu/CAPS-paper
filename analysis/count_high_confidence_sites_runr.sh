#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq


module load R/3.6.0-newgcc

s=caps
log=$HOME/projects/taps/analysis_caps/res/basic_statistics/${s}_count_high_confidence.log

Rscript count_high_confidence.R ../methy_data/${s}_CpG_mods.bed.gz $s 2>&1 | \
tee > $log

s=tapsbeta
s=ps
s=psc

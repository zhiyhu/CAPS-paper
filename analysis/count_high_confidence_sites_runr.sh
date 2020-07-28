#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq
## The following two lines will send an email notification when the job is
## Ended/Aborted/Suspended - Please replace "UserName" with your own username.
#$ -M zyhu
#$ -m aes

module load R/3.6.0-newgcc

s=caps
log=$HOME/projects/taps/analysis_caps/res/basic_statistics/${s}_count_high_confidence.log

Rscript count_high_confidence.R ../methy_data/${s}_CpG_mods.bed.gz $s 2>&1 | \
tee > $log

s=tapsbeta
s=ps
s=psc

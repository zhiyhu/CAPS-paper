#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq
## The following two lines will send an email notification when the job is
## Ended/Aborted/Suspended - Please replace "UserName" with your own username.
#$ -M zyhu
#$ -m aes


## module load bedtools/2.25.0
module load R/3.6.0-newgcc

Rscript 5mc_benchmarking_taps_oxbs.R  2>&1 | \
 tee > logs/5mc_benchmarking_20200726.log
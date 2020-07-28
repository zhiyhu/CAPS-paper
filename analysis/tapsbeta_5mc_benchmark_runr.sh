#!/bin/sh
## The following to run programs in the current working directory
#$ -cwd
## Specify a queue
#$ -q batchq


## module load bedtools/2.25.0
module load R/3.6.0-newgcc

Rscript 5mc_benchmarking_taps_oxbs.R  2>&1 | \
 tee > logs/5mc_benchmarking_20200726.log
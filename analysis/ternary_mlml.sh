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

####################
## prepare tidy data for MLML
Rscript tenary_preprocessing.R

####################
## mlml of TAPS...
dir=/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/mlml
gzip -d ${dir}/tidy_taps_mlml.tsv.gz
gzip -d ${dir}/tidy_tapsbeta_mlml.tsv.gz
gzip -d ${dir}/tidy_caps_mlml.tsv.gz

# mlml -v -u ${dir}/tidy_taps_mlml.tsv -m ${dir}/tidy_tapsbeta_mlml.tsv -h ${dir}/tidy_caps_mlml.tsv | gzip > ${dir}/taps_tapsbeta_caps_mlml_output.tsv.gz
mlml -v                              -m ${dir}/tidy_tapsbeta_mlml.tsv -h ${dir}/tidy_caps_mlml.tsv | gzip > ${dir}/tapsbeta_caps_mlml_output.tsv.gz
# mlml -v -u ${dir}/tidy_taps_mlml.tsv                                  -h ${dir}/tidy_caps_mlml.tsv | gzip > ${dir}/taps_caps_mlml_output.tsv.gz
mlml -v -u ${dir}/tidy_taps_mlml.tsv -m ${dir}/tidy_tapsbeta_mlml.tsv                              | gzip > ${dir}/taps_tapsbeta_mlml_output.tsv.gz

gzip ${dir}/tidy_caps_mlml.tsv
gzip ${dir}/tidy_tapsbeta_mlml.tsv
gzip ${dir}/tidy_taps_mlml.tsv

####################
## mlml of WGBS... 
dir=/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/mlml
gzip -d ${dir}/tidy_wgbs_mlml.tsv.gz
gzip -d ${dir}/tidy_tabseq_mlml.tsv.gz
gzip -d ${dir}/tidy_ace_mlml.tsv.gz

# mlml -v -u ${dir}/tidy_wgbs_mlml.tsv -m ${dir}/tidy_oxbs_mlml.tsv -h ${dir}/tidy_tabseq_mlml.tsv | gzip > ${dir}/wgbs_oxbs_tabseq_mlml_output.tsv.gz
mlml -v -u ${dir}/tidy_wgbs_mlml.tsv -h ${dir}/tidy_tabseq_mlml.tsv | gzip > ${dir}/wgbs_tabseq_mlml_output.tsv.gz
mlml -v -u ${dir}/tidy_wgbs_mlml.tsv -h ${dir}/tidy_ace_mlml.tsv    | gzip > ${dir}/wgbs_ace_mlml_output.tsv.gz

gzip ${dir}/tidy_wgbs_mlml.tsv
gzip ${dir}/tidy_tabseq_mlml.tsv
gzip ${dir}/tidy_ace_mlml.tsv


######################
## tidy up the mlml results
module load R/3.6.0-newgcc
Rscript ternary_filter.R
## filter out conflicts >= 1

####################
## make 1kb windows
module load bedtools/2.25.0 

# genome=/home/obgynae/zyhu/projects/taps/ref/chrom_sizes.tsv
# genomebin=/home/obgynae/zyhu/projects/taps/ref/mm9_1kb_bins.bed.gz
# bedtools makewindows -g $genome -w 1000 | gzip > $genomebin

 # -g <genome>
 #                Genome file size (see notes below).
 #                Windows will be created for each chromosome in the file.
 
 
 # -w <window_size>
 #                Divide each input interval (either a chromosome or a BED interval)
 #                to fixed-sized windows (i.e. same number of nucleotide in each window).
 #                Can be combined with -s <step_size>
 
 
 # (1) The genome file should tab delimited and structured as follows:
 #         <chromName><TAB><chromSize>
 # 
 #        For example, Human (hg19):
 #        chr1    249250621
 #        chr2    243199373
 #        ...
 #        chr18_gl000207_random   4262

####################
## bin mlml results
dir=/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/mlml
genomebin=/home/obgynae/zyhu/projects/taps/ref/mm9_1kb_bins.bed.gz

bed=${dir}/taps_tapsbeta_caps_mlml_output_filtered.tsv.gz
out=${dir}/taps_tapsbeta_caps_mlml_output_filtered_1kbbin.tsv.gz
bedtools map -a $genomebin -b $bed -c 4,5,6 -o mean | gzip > $out

bed=${dir}/tapsbeta_caps_mlml_output_filtered.tsv.gz
out=${dir}/tapsbeta_caps_mlml_output_filtered_1kbbin.tsv.gz
bedtools map -a $genomebin -b $bed -c 4,5,6 -o mean | gzip > $out

bed=${dir}/taps_caps_mlml_output_filtered.tsv.gz
out=${dir}/taps_caps_mlml_output_filtered_1kbbin.tsv.gz
bedtools map -a $genomebin -b $bed -c 4,5,6 -o mean | gzip > $out

bed=${dir}/taps_tapsbeta_mlml_output_filtered.tsv.gz
out=${dir}/taps_tapsbeta_mlml_output_filtered_1kbbin.tsv.gz
bedtools map -a $genomebin -b $bed -c 4,5,6 -o mean | gzip > $out

bed=${dir}/wgbs_ace_mlml_output_filtered.tsv.gz
out=${dir}/wgbs_ace_mlml_output_filtered_1kbbin.tsv.gz
bedtools map -a $genomebin -b $bed -c 4,5,6 -o mean | gzip > $out

bed=${dir}/wgbs_tabseq_mlml_output_filtered.tsv.gz
out=${dir}/wgbs_tabseq_mlml_output_filtered_1kbbin.tsv.gz
bedtools map -a $genomebin -b $bed -c 4,5,6 -o mean | gzip > $out


####################
## intersect mlml 


 
 

 
############ below is obsolete ###########
 
## preprocessing BS-seq file

cd /home/obgynae/zyhu/projects/taps/benchmark/oxbs
zcat GSM4194641_Control_gDNA_100ng.cov.gz | sort -k 1,1 -k2,2n - | gzip > GSM4194641_Control_gDNA_100ng_sorted.cov.gz


mlml -v -u tidy_bs_chrall_mlml2.tsv -h tidy_tabseq_chrall_mlml2.tsv -o bs_tabseq_mlml_output_chrall.tsv
gzip bs_tabseq_mlml_output_chrall.tsv


## TAPS and TAPSbeta
sort -k 1,1 -k2,2n -o tidy_taps_mlml_sorted.tsv tidy_taps_mlml.tsv 


gzip -d tidy_tapsbeta_mlml.tsv.gz
gzip -d tidy_taps_mlml.tsv.gz

mlml -v -u tidy_taps_mlml.tsv -m tidy_tapsbeta_mlml.tsv -o taps_tapsbeta_mlml_output.tsv
gzip taps_tapsbeta_mlml_output.tsv
gzip tidy_tapsbeta_mlml.tsv
gzip tidy_taps_mlml.tsv



## Zhiyuan Hu
## 4 july 2020
## Using mlml from the MethPipe
## mlml parameters:
##  -o, -output     Name of output file (default: stdout) 
##  -u, -bsseq      Name of input BS-Seq methcounts file 
##  -h, -tabseq     Name of input TAB-Seq methcounts file 
##  -m, -oxbsseq    Name of input oxBS-Seq methcounts file 
##  -t, -tolerance  EM convergence threshold. Default 1e-10 
##  -a, -alpha      significance level of binomial test for each site. Default 
##                  0.05 
##  -v, -verbose    print run statistics 

## Input data format
## chr1    1       +       CpG     0.363636363636  11
## chr1    3       +       CpG     1.0     7
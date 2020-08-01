## Zhiyuan Hu
## last modified 1 Aug 2020

module load bedtools/2.25.0 

beda=$HOME/projects/taps/analysis_caps/tidy_data/5hmC_caps_mlml_ace_tab_mods.unfiltered.bed.gz
beda2=$HOME/projects/taps/analysis_caps/tidy_data/5hmC_caps_mlml_ace_tab_mods.unfiltered_tidy.bed.gz
zcat $beda | grep -v 'strand' | gzip > $beda2 # remove the header

for region in bins left_4kb_flank_20wins right_4kb_flank_20wins
do
dir=$HOME/projects/taps/analysis_caps/tidy_data
beda=${dir}/5hmC_caps_mlml_ace_tab_mods.unfiltered_tidy.bed.gz
bedb=$HOME/projects/taps/ref/ucsc/cpgIslandExt_${region}.*.gz

for i in 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 X Y
do
out=${dir}/cpg_island/5hmC_caps_mlml_ace_tab_cpgIsland_${region}_annotated_chr${i}.bed.gz
zcat $bedb | grep 'chr'$i | bedtools intersect -wo -a $beda -b - | gzip > $out
echo $i
done
done

####### run R

module load R/3.6.0-newgcc

Rscript CpG_island_analysis.R



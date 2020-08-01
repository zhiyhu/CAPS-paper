#!/usr/bin/env Rscript
# preprocessing data for the IGV
# Zhiyuan Hu
# 6 july 2020

# Output file: tidy_data/igv/{method}_modlevel.bed.gz

library(data.table)
suppressPackageStartupMessages(library(parallel))
library(MASS)
suppressPackageStartupMessages(library(dplyr))

fn_taps <- "$HOME/projects/taps/analysis_caps/methy_data/taps_CpG_mods.bed.gz"
fn_beta <- "$HOME/projects/taps/analysis_caps/methy_data/tapsbeta_CpG_mods.bed.gz"
fn_caps <- "$HOME/projects/taps/analysis_caps/methy_data/caps_CpG_mods.bed.gz"

fn_mod <- c("$HOME/projects/taps/analysis_caps/tidy_data/igv/taps_modlevel.bed.gz",
            "$HOME/projects/taps/analysis_caps/tidy_data/igv/tapsbeta_modlevel.bed.gz",
            "$HOME/projects/taps/analysis_caps/tidy_data/igv/caps_modlevel.bed.gz")

## the positions need to be match

taps <- fread(fn_taps, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(taps) <- c("chr","start","end","mod_level","mod","unmod")
taps <- taps[taps$chr %in% paste("chr", 1:19, sep = ""),]
taps <- taps[taps$mod + taps$unmod > 0,]

beta <- fread(fn_beta, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(beta) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")
beta <- beta[beta$chr %in% paste("chr", 1:19, sep = ""),]
beta <- beta[beta$snv == "No",]
beta <- beta[beta$mod + beta$unmod > 0,]

caps <- fread(fn_caps, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(caps) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")
caps <- caps[caps$chr %in% paste("chr", 1:19, sep = ""),]
caps <- caps[caps$snv == "No",]
caps <- caps[caps$mod + caps$unmod > 0,]

print("tidy up")

tidy_beta <- data.frame(chr     = beta$chr,
                        start   = beta$start,
                        end     = beta$end,
                        mod_level = beta$mod/(beta$mod + beta$unmod))

tidy_caps <- data.frame(chr     = caps$chr,
                        start   = caps$start,
                        end     = caps$end,
                        mod_level = caps$mod/(caps$mod + caps$unmod))

tidy_taps <- data.frame(chr     = taps$chr,
                        start   = taps$start,
                        end     = taps$end,
                        mod_level = taps$mod/(taps$mod + taps$unmod))


fwrite(x = tidy_taps[,c(1:4)], fn_mod[1], sep = "\t", col.names = F, row.names = F)
fwrite(x = tidy_beta[,c(1:4)], fn_mod[2], sep = "\t", col.names = F, row.names = F)
fwrite(x = tidy_caps[,c(1:4)], fn_mod[3], sep = "\t", col.names = F, row.names = F)


## get TAB-seq and ACE-seq dat

fn <- "$HOME/projects/taps/analysis_caps/tidy_data/5hmC_caps_mlml_ace_tab_mods.unfiltered.bed.gz"
fn_mod <- c("$HOME/projects/taps/analysis_caps/tidy_data/igv/tabseq_modlevel.bed.gz",
            "$HOME/projects/taps/analysis_caps/tidy_data/igv/ace_modlevel.bed.gz")

mod <- fread(fn, stringsAsFactors = F, data.table = F, header = T, sep = "\t", verbose = T, nThread = 4)

## tabseq
tab <- mod[,c("chr","start","end","tab.mod","tab.unmod")]
tab$depth <- tab$tab.mod + tab$tab.unmod
tab <- tab[tab$depth > 0,]
tab$modlevel <- tab$tab.mod/tab$depth

fwrite(x = tab[,c("chr","start","end","modlevel")], fn_mod[1], sep = "\t", col.names = F, row.names = F)

## ace
ace <- mod[,c("chr","start","end","ace.mod","ace.unmod")]
ace$depth <- ace$ace.mod + ace$ace.unmod
ace <- ace[ace$depth > 0,]
ace$modlevel <- ace$ace.mod/ace$depth

fwrite(x = ace[,c("chr","start","end","modlevel")], fn_mod[2], sep = "\t", col.names = F, row.names = F)


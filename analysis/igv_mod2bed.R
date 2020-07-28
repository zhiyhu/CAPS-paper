#!/usr/bin/env Rscript
# preprocessing data for the IGV
# Zhiyuan Hu
# 6 july 2020

library(data.table)
suppressPackageStartupMessages(library(parallel))
library(MASS)
suppressPackageStartupMessages(library(dplyr))

fn_taps <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/taps_CpG_mods.bed.gz"
fn_beta <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/tapsbeta_CpG_mods.bed.gz"
fn_caps <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/caps_CpG_mods.bed.gz"

fn_mod <- c("/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/taps_modlevel.bed.gz",
            "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/tapsbeta_modlevel.bed.gz",
            "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/caps_modlevel.bed.gz")
fn_dep <- c("/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/taps_depth.bed.gz",
            "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/tapsbeta_depth.bed.gz",
            "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/caps_depth.bed.gz")
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
                        mod_level = beta$mod/(beta$mod + beta$unmod),
                        depth     = beta$mod + beta$unmod)

tidy_caps <- data.frame(chr     = caps$chr,
                        start   = caps$start,
                        end     = caps$end,
                        mod_level = caps$mod/(caps$mod + caps$unmod),
                        depth     = caps$mod + caps$unmod)

tidy_taps <- data.frame(chr     = taps$chr,
                        start   = taps$start,
                        end     = taps$end,
                        mod_level = taps$mod/(taps$mod + taps$unmod),
                        depth     = taps$mod + taps$unmod)


fwrite(x = tidy_taps[,c(1:4)], fn_mod[1], sep = "\t", col.names = F, row.names = F)
fwrite(x = tidy_beta[,c(1:4)], fn_mod[2], sep = "\t", col.names = F, row.names = F)
fwrite(x = tidy_caps[,c(1:4)], fn_mod[3], sep = "\t", col.names = F, row.names = F)

fwrite(x = tidy_taps[,c(1:3,5)], fn_dep[1], sep = "\t", col.names = F, row.names = F)
fwrite(x = tidy_beta[,c(1:3,5)], fn_dep[2], sep = "\t", col.names = F, row.names = F)
fwrite(x = tidy_caps[,c(1:3,5)], fn_dep[3], sep = "\t", col.names = F, row.names = F)

## get TAB-seq and ACE-seq dat

fn <- "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/5hmC_caps_mlml_ace_tab_mods.unfiltered.bed.gz"
fn_mod <- c("/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/tabseq_modlevel.bed.gz",
            "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/tabseq_depth.bed.gz",
            "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/ace_modlevel.bed.gz",
            "/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/ace_depth.bed.gz")

mod <- fread(fn, stringsAsFactors = F, data.table = F, header = T, sep = "\t", verbose = T, nThread = 4)

## tabseq
tab <- mod[,c("chr","start","end","tab.mod","tab.unmod")]
tab$depth <- tab$tab.mod + tab$tab.unmod
tab <- tab[tab$depth > 0,]
tab$modlevel <- tab$tab.mod/tab$depth

fwrite(x = tab[,c("chr","start","end","modlevel")], fn_mod[1], sep = "\t", col.names = F, row.names = F)
fwrite(x = tab[,c("chr","start","end","depth")], fn_mod[2], sep = "\t", col.names = F, row.names = F)


## ace
ace <- mod[,c("chr","start","end","ace.mod","ace.unmod")]
ace$depth <- ace$ace.mod + ace$ace.unmod
ace <- ace[ace$depth > 0,]
ace$modlevel <- ace$ace.mod/ace$depth

fwrite(x = ace[,c("chr","start","end","modlevel")], fn_mod[3], sep = "\t", col.names = F, row.names = F)
fwrite(x = ace[,c("chr","start","end","depth")], fn_mod[4], sep = "\t", col.names = F, row.names = F)

##############################################################
## extact positions where three methods captured high methylation

mod$caps.depth <- mod$caps.unmod + mod$caps.mod
mod$ace.depth <- mod$ace.unmod + mod$ace.mod
mod$tab.depth <- mod$tab.unmod + mod$tab.mod

mod <- mod[mod$caps.depth >=3 &
             mod$ace.depth >= 3 &
             mod$tab.depth >= 3, ]

mod$caps.modlevel <- mod$caps.mod / mod$caps.depth
mod$ace.modlevel  <- mod$ace.mod  / mod$ace.depth
mod$tab.modlevel  <- mod$tab.mod  / mod$tab.depth


sum(mod$caps.modlevel > 0.2 &
      mod$ace.modlevel> 0.2 &
      mod$tab.modlevel> 0.2)
# 5157

sum(mod$caps.modlevel > 0.3 &
      mod$ace.modlevel> 0.3 &
      mod$tab.modlevel> 0.3)
# 322

sum(mod$caps.modlevel > 0.4 &
      mod$ace.modlevel> 0.4 &
      mod$tab.modlevel> 0.4)
# 9

fwrite(mod[mod$caps.modlevel > 0.4 &
             mod$ace.modlevel> 0.4 &
             mod$tab.modlevel> 0.4,1:4],"/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/5hmc_caps_ace_tab_sites_0.4.bed.gz", sep = "\t", col.names = T, row.names = F)


fwrite(mod[mod$caps.modlevel > 0.3 &
             mod$ace.modlevel> 0.3 &
             mod$tab.modlevel> 0.3,1:4],"/home/obgynae/zyhu/projects/taps/analysis_caps/tidy_data/igv/5hmc_caps_ace_tab_sites_0.3.bed.gz", sep = "\t", col.names = T, row.names = F)




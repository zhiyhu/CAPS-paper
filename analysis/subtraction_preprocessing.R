#!/usr/bin/env Rscript
# preprocessing data for the ternary plot
# Zhiyuan Hu
# 4 july 2020
# last modified 18 Nov 2020

library(data.table)
suppressPackageStartupMessages(library(parallel))
library(MASS)
suppressPackageStartupMessages(library(dplyr))

## set up ----
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
} else {
  os <- "No"
}

#---------------------------#
#  TAPS, TAPSbeta, CAPS     #
#---------------------------#

fn_taps <- "$HOME/projects/taps/analysis_caps/methy_data/taps_CpG_mods.bed.gz"
fn_beta <- "$HOME/projects/taps/analysis_caps/methy_data/tapsbeta_CpG_mods.bed.gz"
fn_caps <- "$HOME/projects/taps/analysis_caps/methy_data/caps_CpG_mods.bed.gz"
chrs <- paste("chr", c(1:19, "X","Y"), sep = "")
fn_out <- c("$HOME/projects/taps/analysis_caps/tidy_data/mlml/tidy_taps_mlml.tsv.gz",
              "$HOME/projects/taps/analysis_caps/tidy_data/mlml/tidy_tapsbeta_mlml.tsv.gz",
              "$HOME/projects/taps/analysis_caps/tidy_data/mlml/tidy_caps_mlml.tsv.gz")

## the positions need to be match

taps <- fread(fn_taps, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(taps) <- c("chr","start","end","mod_level","mod","unmod")

beta <- fread(fn_beta, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(beta) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")
beta <- beta[beta$snv == "No",]

caps <- fread(fn_caps, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(caps) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")
caps <- caps[caps$snv == "No",]

tidy_beta <- data.frame(chr     = beta$chr,
                        pos     = beta$start,
                        strand  = "+",
                        context = "CpG",
                        mod_level = beta$mod_level,
                        depth     = beta$mod + beta$unmod)

tidy_caps <- data.frame(chr     = caps$chr,
                        pos     = caps$start,
                        strand  = "+",
                        context = "CpG",
                        mod_level = caps$mod_level,
                        depth     = caps$mod + caps$unmod)

tidy_taps <- data.frame(chr     = taps$chr,
                        pos     = taps$start,
                        strand  = "+",
                        context = "CpG",
                        mod_level = taps$mod_level,
                        depth     = taps$mod + taps$unmod)

tidy_taps <- tidy_taps[tidy_taps$depth >= 5,]
tidy_beta <- tidy_beta[tidy_beta$depth >= 5,]
tidy_caps <- tidy_caps[tidy_caps$depth >= 5,]

tidy_data <- inner_join(tidy_taps[,1:2], tidy_beta[,1:2], by = c("chr","pos"))
tidy_data <- inner_join(tidy_data[,1:2], tidy_caps[,1:2], by = c("chr","pos"))

tidy_beta <- inner_join(tidy_beta, tidy_data, by = c("chr","pos"))
tidy_caps <- inner_join(tidy_caps, tidy_data, by = c("chr","pos"))
tidy_taps <- inner_join(tidy_taps, tidy_data, by = c("chr","pos"))

tidy_taps <- tidy_taps[order(tidy_taps$chr, tidy_taps$pos),]
tidy_beta <- tidy_beta[order(tidy_beta$chr, tidy_beta$pos),]
tidy_caps <- tidy_caps[order(tidy_caps$chr, tidy_caps$pos),]

if(os != "Darwin") {
  fwrite(x = tidy_taps, fn_out[1], sep = "\t", col.names = F)
  fwrite(x = tidy_beta, fn_out[2], sep = "\t", col.names = F)
  fwrite(x = tidy_caps, fn_out[3], sep = "\t", col.names = F)
}


## BS, ACE-seq and TAB-seq -----

fn_wgbs <- "$HOME/projects/taps/analysis_caps/methy_data/wgbs_CpG_mods.bed.gz"
fn_ace <- "$HOME/projects/taps/analysis_caps/methy_data/ace_CpG_mods_noheader.bed.gz"
fn_tab <- "$HOME/projects/taps/analysis_caps/methy_data/tabseq_CpG_mods.bed.gz"
chrs <- paste("chr", c(1:19, "X","Y"), sep = "")
fn_out <- c("$HOME/projects/taps/analysis_caps/tidy_data/mlml/tidy_wgbs_mlml.tsv.gz",
            "$HOME/projects/taps/analysis_caps/tidy_data/mlml/tidy_ace_mlml.tsv.gz",
            "$HOME/projects/taps/analysis_caps/tidy_data/mlml/tidy_tabseq_mlml.tsv.gz")

## the positions need to be match

wgbs <- fread(fn_wgbs, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(wgbs) <- c("chr","start","end","mod_level","mod","unmod")

ace <- fread(fn_ace, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(ace) <- c("chr","start","end","mod_level","depth","strand")

tab <- fread(fn_tab, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(tab) <- c("chr","start","end","mod_level","mod","unmod")

tidy_ace <- data.frame(chr     = ace$chr,
                       pos     = ace$start,
                       strand  = "+",
                       context = "CpG",
                       mod_level = ace$mod_level,
                       depth     = ace$depth)

tidy_tab <- data.frame(chr     = tab$chr,
                        pos     = tab$start,
                        strand  = "+",
                        context = "CpG",
                        mod_level = tab$mod/(tab$mod + tab$unmod),
                        depth     = tab$mod + tab$unmod)

tidy_wgbs <- data.frame(chr     = wgbs$chr,
                        pos     = wgbs$start,
                        strand  = "+",
                        context = "CpG",
                        mod_level = wgbs$mod/(wgbs$mod + wgbs$unmod),
                        depth     = wgbs$mod + wgbs$unmod)

tidy_wgbs <- tidy_wgbs[tidy_wgbs$depth >= 5,]
tidy_ace <- tidy_ace[   tidy_ace$depth >= 5,]
tidy_tab <- tidy_tab[   tidy_tab$depth >= 5,]

tidy_data <- inner_join(tidy_wgbs[,1:2], tidy_ace[,1:2], by = c("chr","pos"))
tidy_data <- inner_join(tidy_data[,1:2], tidy_tab[,1:2], by = c("chr","pos"))


tidy_wgbs <- inner_join(tidy_wgbs, tidy_data, by = c("chr","pos"))
tidy_ace <- inner_join(tidy_ace, tidy_data, by = c("chr","pos"))
tidy_tab <- inner_join(tidy_tab, tidy_data, by = c("chr","pos"))

tidy_wgbs <- tidy_wgbs[order(tidy_wgbs$chr, tidy_wgbs$pos),]
tidy_ace <- tidy_ace[order(tidy_ace$chr, tidy_ace$pos),]
tidy_tab  <- tidy_tab [order(tidy_tab$chr, tidy_tab$pos),]


fwrite(x = tidy_wgbs, fn_out[1], sep = "\t", col.names = F)
fwrite(x = tidy_ace, fn_out[2], sep = "\t", col.names = F)
fwrite(x = tidy_tab, fn_out[3], sep = "\t", col.names = F)

#!/usr/bin/env Rscript
## preprocessing TAPS TAPSbeta for mlml
## Zhiyuan Hu
## 11 May 2020
## last updated 26 july 2020

library(data.table)

## set up ----
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
} else {
  os <- "No"
}

fn_taps <- "$HOME/projects/taps/analysis_caps/methy_data/taps_CpG_mods.bed.gz"
fn_beta <- "$HOME/projects/taps/analysis_caps/methy_data/tapsbeta_CpG_mods.bed.gz"
chrs <- paste("chr", c(1:19, "X","Y"), sep = "")
fn_out <- c("$HOME/projects/taps/caps/R/tidy_data/mlml/tidy_taps_mlml.tsv.gz",
              "$HOME/projects/taps/caps/R/tidy_data/mlml/tidy_tapsbeta_mlml.tsv.gz")


## the positions need to be match

taps <- fread(fn_taps, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4, )
colnames(taps) <- c("chr","start","end","mod_level","mod","unmod")

beta <- fread(fn_beta, stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
colnames(beta) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")

tidy_beta <- data.frame(chr     = beta$chr,
                        pos     = beta$start,
                        strand  = "+",
                        context = "CpG",
                        mod_level = beta$mod_level,
                        depth     = beta$mod + beta$unmod)

tidy_taps <- data.frame(chr     = taps$chr,
                        pos     = taps$start,
                        strand  = "+",
                        context = "CpG",
                        mod_level = taps$mod_level,
                        depth     = taps$mod + taps$unmod)

tidy_taps <- tidy_taps[tidy_taps$depth >= 3,]
tidy_beta <- tidy_beta[tidy_beta$depth >= 3,]

tidy_data <- merge(tidy_taps[,1:2], tidy_beta[,1:2], by = c("chr","pos"))

tidy_beta <- merge(tidy_beta, tidy_data, by = c("chr","pos"))
tidy_taps <- merge(tidy_taps, tidy_data, by = c("chr","pos"))


if(os != "Darwin") {
  fwrite(x = tidy_taps, fn_out[1], sep = "\t")
  fwrite(x = tidy_beta, fn_out[2], sep = "\t")
}
#!/usr/bin/env Rscript
## colocolization analysis with genomic regulatory elements of 5hmC
## Zhiyuan Hu
## 20 May 2020
## last modified: 28 Jun 2020
## ran remotely

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

depth_cutoff <- 5
p_cutoff <- 0.05
ncores <- 4

## parsing arguments
options <- commandArgs(trailingOnly = TRUE)

## false positive rate
fp_rate <- read.table("false_positive_table.txt", header = T, stringsAsFactors = F)
fp_rate$methods <- toupper(fp_rate$methods)

## part 2 ----

fn <- options[1]
method <- options[2]
method <- toupper(method)

print(paste0("fn is ", fn))
print(paste0("method is ", method))

caps_all <- fread(fn, header = F, nThread = 4, stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
colnames(caps_all) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")

## filter variants and depth
caps_all <- caps_all[caps_all$snv == "No",] ## filter by SNV
caps_all$depth <- caps_all$mod + caps_all$unmod
caps_all <- caps_all[caps_all$depth >= depth_cutoff,] # filter by depth
caps_all <- caps_all[,c("chr","start","end","mod_level","mod","unmod","depth")]

## binomial test ----
binom_t <- function(depth, mod, p) { # false positive rate/probability
  p_value <- mcmapply(as.numeric(mod), as.numeric(depth), 
                      FUN = function(x, y) {
                        if(!is.na(x) & !is.na(y)) {
                          if(y > 0) {
                            return(binom.test(x = x, n = y, p = p, alternative = "greater")$p.value)
                          } else {
                            return(1)
                          }
                        } else {
                          return(1)
                        }
                        
                      }, 
                      mc.cores = ncores)
  return(p_value)
}

# print("binomial test caps")

fp <- fp_rate$conversion[fp_rate$methods == method]
p <- binom_t(depth = caps_all$depth, mod = caps_all$mod, p = fp)
rm(caps_all)

fdr <- p.adjust(p, method = "BH")
n <- sum(fdr <= p_cutoff)

print(paste0("before test n = ", length(p)))
print(paste0("after test n = ", n))



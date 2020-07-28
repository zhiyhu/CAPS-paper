#!/usr/bin/env Rscript
## Benchmark CAPS and ACE, TAB-seq
## author: Zhiyuan Hu
## created 25 Apr 2020
## last modified 25 Apr 2020

library(data.table)
library(parallel)
library(MASS)
library(ggplot2)
library(ggplotify)

local = T

# date  = format(Sys.time(), "%d%b%Y") # date
# output  = "" # output file name
# assays = c("caps", "ace", "tabseq") # assays  
ncores = 4 # n CPUs

## CAPS, ACE, TABseq ----
## read in bedtools map output

df_mod <- read.delim("../../data/tidy_data/5hmC_caps_mlml_ace_tab_mods_10kb_bin.bed.gz", as.is = T, header = F)
colnames(df_mod) <- c("chr","start","end",'caps_mod','caps_unmod','ace_mod','ace_unmod','tab_mod','tab_unmod','mlml_mod','mlml_unmod',"caps_mod_raw", "ace_mod_raw","tab_mod_raw")

dim(df_mod)
# [1] 265500     12
# filter out uncovered bins
df_mod <- df_mod[ !  (is.na(df_mod$caps_mod)   & 
                     is.na(df_mod$caps_unmod) &
                     is.na(df_mod$ace_mod)    &
                     is.na(df_mod$ace_unmod)  &
                     is.na(df_mod$tab_mod)    &
                     is.na(df_mod$tab_unmod)  &
                     is.na(df_mod$mlml_mod)   &
                     is.na(df_mod$mlml_unmod)),]
dim(df_mod)
# [1] 254357     14


for(i in 4:ncol(df_mod)) {
  df_mod[,i] <- as.numeric(df_mod[,i])
}

# calculate raw signals
df_mod$caps_signal <- df_mod$caps_mod/(df_mod$caps_mod + df_mod$caps_unmod)
df_mod$ace_signal  <- df_mod$ace_mod /(df_mod$ace_mod  + df_mod$ace_unmod )
df_mod$tab_signal  <- df_mod$tab_mod /(df_mod$tab_mod  + df_mod$tab_unmod )
df_mod$mlml_signal <- df_mod$mlml_mod/(df_mod$mlml_mod + df_mod$mlml_unmod)

df_mod$raw_caps_signal <- df_mod$caps_mod_raw/(df_mod$caps_mod_raw  + df_mod$caps_unmod)
df_mod$raw_ace_signal  <- df_mod$ace_mod_raw /( df_mod$ace_mod_raw  + df_mod$ace_unmod )
df_mod$raw_tab_signal  <- df_mod$tab_mod_raw /( df_mod$tab_mod_raw  + df_mod$tab_unmod )

# choose different cutoffs -----  

i_cutoff <-  500

hist((df_mod$caps_mod + df_mod$caps_unmod), 200) ## 200 is a turning point for CAPS

hist((df_mod$caps_mod + df_mod$caps_unmod), 500, xlim = c(0, 2000)) ## 200 is a turning point for CAPS

df_mod2 <- df_mod[(df_mod$caps_mod + df_mod$caps_unmod) >= i_cutoff & (df_mod$caps_mod + df_mod$caps_unmod) < 5000 &
                    (df_mod$ace_mod  + df_mod$ace_unmod) >= i_cutoff & (df_mod$ace_mod  + df_mod$ace_unmod) < 5000 &
                    (df_mod$tab_mod  + df_mod$tab_unmod) >= i_cutoff & (df_mod$tab_mod  + df_mod$tab_unmod) < 5000, ]


dim(df_mod)
dim(df_mod2)
# [1] 209681     21

df_compare <- data.frame(assay1 = c("CAPS","CAPS","TAB-seq","mlml","mlml"),
                         assay2 = c("ACE-seq","TAB-seq", "ACE-seq","ACE-seq","TAB-seq"),
                         person_cor = NA,
                         spearman_cor = NA,
                         n_bins = nrow(df_mod2),
                         depth_cutoff_for_bins = i_cutoff)

idx <- !is.na(df_mod2$raw_caps_signal) & !is.na(df_mod2$raw_ace_signal)
df_compare$person_cor[1] <- cor(df_mod2$raw_caps_signal[idx],  df_mod2$raw_ace_signal[idx], method = "pearson")
df_compare$spearman_cor[1] <- cor(df_mod2$raw_caps_signal[idx],  df_mod2$raw_ace_signal[idx], method = "spearman")

p1 <- as.ggplot(~smoothScatter(x = df_mod2$raw_caps_signal[idx],  y = df_mod2$raw_ace_signal[idx], xlab = "CAPS", ylab = "ACE-seq",
                               xlim = c(0, 0.2), ylim = c(0,0.2))) 

idx <- !is.na(df_mod2$raw_caps_signal) & !is.na(df_mod2$raw_tab_signal)
df_compare$person_cor[2] <-  cor(df_mod2$raw_caps_signal[idx], df_mod2$raw_tab_signal[idx])
df_compare$spearman_cor[2] <-  cor(df_mod2$raw_caps_signal[idx], df_mod2$raw_tab_signal[idx], method = "spearman")
p2 <- as.ggplot(~smoothScatter(x = df_mod2$raw_caps_signal[idx],  y = df_mod2$raw_tab_signal[idx], xlab = "CAPS", ylab = "TAB-seq",
                               xlim = c(0, 0.2), ylim = c(0,0.2))) 


idx <- !is.na(df_mod2$raw_ace_signal) & !is.na(df_mod2$raw_tab_signal)
df_compare$person_cor[3] <-  cor(df_mod2$raw_ace_signal[idx],  df_mod2$raw_tab_signal[idx])
df_compare$spearman_cor[3] <-  cor(df_mod2$raw_ace_signal[idx],  df_mod2$raw_tab_signal[idx], method = "spearman")
p3 <- as.ggplot(~smoothScatter(x = df_mod2$raw_ace_signal[idx], y = df_mod2$raw_tab_signal[idx], xlab = "ACE-seq", ylab = "TAB-seq",
                               xlim = c(0, 0.2), ylim = c(0,0.2))) 


idx <- !is.na(df_mod2$mlml_signal) & !is.na(df_mod2$ace_signal)
df_compare$n_bins[4] <- sum(idx)
df_compare$person_cor[4] <- cor(df_mod2$mlml_signal[idx],  df_mod2$raw_ace_signal[idx])
df_compare$spearman_cor[4] <- cor(df_mod2$mlml_signal[idx],  df_mod2$raw_ace_signal[idx], method = "spearman")
p4 <- as.ggplot(~smoothScatter(x = df_mod2$mlml_signal[idx],  y = df_mod2$raw_ace_signal[idx], xlab = "MLML", ylab = "ACE-seq",
                xlim = c(0, 0.2), ylim = c(0,0.2)))

idx <- !is.na(df_mod2$mlml_signal) & !is.na(df_mod2$tab_signal)
df_compare$n_bins[5] <- sum(idx)
df_compare$person_cor[5]   <- cor(df_mod2$mlml_signal[idx], df_mod2$raw_tab_signal[idx])
df_compare$spearman_cor[5] <- cor(df_mod2$mlml_signal[idx], df_mod2$raw_tab_signal[idx], method = "spearman")
p5 <- as.ggplot(~smoothScatter(x = df_mod2$mlml_signal[idx],  y = df_mod2$raw_tab_signal[idx], xlab = "MLML", ylab = "TAB-seq",
                               xlim = c(0, 0.2), ylim = c(0,0.2)))

# idx <- !is.na(df_mod2$mlml_signal) & !is.na(df_mod2$caps_signal)
# df_compare$n_bins[5] <- sum(idx)
# df_compare$person_cor[5]   <- cor(df_mod2$mlml_signal[idx], df_mod2$raw_tab_signal[idx])
# df_compare$spearman_cor[5] <- cor(df_mod2$mlml_signal[idx], df_mod2$raw_tab_signal[idx], method = "spearman")
# p6 <- as.ggplot(~smoothScatter(x = df_mod2$mlml_signal[idx],  y = df_mod2$raw_caps_signal[idx], xlab = "MLML", ylab = "CAPS"))

cowplot::plot_grid(p1,p2,p3,p4,p5, ncol = 3)
ggsave("plots/benchmark_caps200515/caps_mlml_ace_tabseq.smoothScatter_adjusted.pdf", width = 8.3, height = 6)
# write.csv(df_compare, "plots/benchmark_caps200515/correlation_results.csv", row.names = F)


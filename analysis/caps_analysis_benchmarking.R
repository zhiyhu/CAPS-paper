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





### ### ### ### ### ### ### ###
########### below is obsolete
### ### ### ### ### ### ### ###

## there is a huge difference between caps and ace in some bins regarding coverage
pdf("plots/benchmark_caps200426/diff_in_bins_coverage.pdf")
ggplot(df_mod2, aes(caps_mod + caps_unmod, ace_mod + ace_unmod)) + 
  geom_point(pch = 21, alpha = 0.5) + 
  geom_abline(intercept = 1000, slope = 0.9)+ 
  geom_abline(intercept = -2000, slope = 0.9)
dev.off()

## write the bins with difference
write.table(x = df_mod2[which(df_mod2$caps_mod+ df_mod2$caps_unmod  > df_mod2$ace_mod+ df_mod2$ace_unmod + 3000),], 
            file = "plots/benchmark_caps200426/bins_with_more_caps_coverage.tsv",
            quote = F, sep = "\t", row.names = F)


## filter out this different bins


summary(lm(df_mod2$caps_cov ~ df_mod2$ace_cov))
# intercept 2.443e+02
# slope 1.132e+00

df_mod2$outliers <- F
df_mod2$outliers[df_mod2$caps_cov/df_mod2$ace_cov >= 2 | df_mod2$ace_cov/df_mod2$caps_cov >= 2] <- T

ggplot(df_mod2, aes(caps_cov, ace_cov, color = outliers)) + 
  geom_point(pch = 21, alpha = 0.5) + 
  geom_abline(intercept = 0, slope =2)+ 
  geom_abline(intercept = 0, slope = 1/2)

df_mod2 <- df_mod2[df_mod2$outliers == F,]

## look into the difference via bed_intersect.sh


### bin manually

# d <- MASS::kde2d(x = df_mod2$caps_signal, y = df_mod2$ace_signal, n = 299, lims = c(0, 0.3, 0, 0.3))
# mat_d <- data.frame(d$x, d$z)
# colnames(mat_d) <- c("x", d$y)
# mat_d <- tidyr::gather(mat_d, y, level, "0":"0.3", factor_key=TRUE) # http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/

# df_mod2$caps_sig_bin <- cut(x = df_mod2$caps_signal, breaks = seq(0, 0.2, 0.001))
# df_mod2$ace_sig_bin  <- cut(x = df_mod2$ace_signal,  breaks = seq(0, 0.2, 0.001))
# df_mod2$tab_sig_bin  <- cut(x = df_mod2$tab_signal,  breaks = seq(0, 0.2, 0.001))
# 
# n_bin <- length(levels(df_mod2$caps_sig_bin))
# y_bin <- data.frame(id = paste("bin", 1:n_bin, sep = ""),
#                     range = levels(df_mod2$ace_sig_bin))
# y_bin$low <- sapply(as.character(y_bin$range), function(x) return(unlist(strsplit(x, "[,]"))[1]))
# y_bin$low <- as.numeric(gsub(y_bin$low, pattern = "[(]", replacement = ""))
# y_bin$upper <- sapply(as.character(y_bin$range), function(x) return(unlist(strsplit(x, "[,]"))[2]))
# y_bin$upper  <- as.numeric(gsub(y_bin$upper , pattern = "]", replacement = ""))
# y_bin$mid <- (y_bin$low + y_bin$upper)/2
# 
# x_bin <- data.frame(id = paste("bin", 1:n_bin, sep = ""),
#                     range = levels(df_mod2$caps_sig_bin))
# x_bin$low <- sapply(as.character(x_bin$range), function(x) return(unlist(strsplit(x, "[,]"))[1]))
# x_bin$low <- as.numeric(gsub(x_bin$low, pattern = "[(]", replacement = ""))
# x_bin$upper <- sapply(as.character(x_bin$range), function(x) return(unlist(strsplit(x, "[,]"))[2]))
# x_bin$upper  <- as.numeric(gsub(x_bin$upper , pattern = "]", replacement = ""))
# x_bin$mid <- (x_bin$low + x_bin$upper)/2
# 
# den <- table(df_mod2$caps_sig_bin, df_mod2$ace_sig_bin)
# den <- data.frame(cbind(x_bin$mid, den))
# colnames(den) <- c("x", as.character(y_bin$mid))
# mat_d <- tidyr::gather(den, y, level, 2:ncol(den), factor_key = F) # http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
# 
# ## gg raster
# ggplot(mat_d[mat_d$level > 0,], aes(x, as.numeric(as.character(y)))) +  theme_linedraw() +
#   geom_tile(aes(fill = log2(level + 1)) ) + scale_fill_viridis(option = "C") + 
#   xlab("CAPS raw signal") + 
#   ylab("ACE-seq raw signal") +
#   xlim(0, max(df_mod2$caps_signal) + 0.002) +
#   ylim(0,  max(df_mod2$caps_signal) + 0.002)
# ggsave(paste("plots/benchmark_caps200426/raster_10kbBin_rawSig_CAPS_ACE_", 
#              nrow(df_mod2), "bins_cor", round(cor(df_mod2$caps_signal, df_mod2$ace_signal), digits = 3), 
#              "_cutoff", i_cutoff, 
#              ".pdf", sep = ""))




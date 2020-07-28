#!/usr/bin/env Rscript
## Tenary diagram 
## author: Zhiyuan Hu
## created 23 Apr 2020
## last modified 4 july 2020
## creat ternary plot


### ternary plot


library('Ternary')
library(RColorBrewer)

TernaryDiagramPrepare <- function(filename) {
  taps <- fread(filename, 
                stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
  colnames(taps) <- c('chr','start','end','mc','hmc','c')
  taps$c <- as.numeric(taps$c)
  taps <- taps[!is.na(taps$c),]
  taps$mc <- as.numeric(taps$mc)
  taps$hmc <- as.numeric(taps$hmc)
  # taps <- taps[taps$c >0 & taps$hmc >0  & taps$mc > 0, ]
  m <- as.matrix(taps[,c(5,6,4)]) # in order of hmc, c, mc
  return(m)
}

colfunc <- colorRampPalette(colors = c("white",brewer.pal(9,"BuPu")[3:9], "black"))
pdf("plots/ternary_plot20200705/color_gradient.pdf")
plot(rep(1,100),col=(colfunc(100)), pch=15,cex=2)
dev.off()

col <- alpha(colfunc(500),0.8)
# colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))

### TAPS
# 
# taps <- TernaryDiagramPrepare(filename = "../../data/mlml_output/taps_tapsbeta_caps_mlml_output_filtered_1kbbin.tsv.gz")
# density <- TernaryDensity(taps, resolution = 25L)
# pdf("plots/ternary_plot20200705/taps_tapsbeta_caps_mlml.pdf")
# TernaryPlot(axis.labels = seq(0, 100, by = 10), atip = "5hmC",btip = "C",ctip = "5mC")
# ColourTernary(density, spectrum = col)
# dev.off()

# taps <- TernaryDiagramPrepare(filename = "../../data/mlml_output/taps_caps_mlml_output_filtered_1kbbin.tsv.gz")
# density <- TernaryDensity(taps, resolution = 25L)
# pdf("plots/ternary_plot20200705/taps_caps_mlml.pdf")
# TernaryPlot(axis.labels = seq(0, 100, by = 10), atip = "5hmC",btip = "C",ctip = "5mC")
# ColourTernary(density, spectrum = col)
# dev.off()

taps <- TernaryDiagramPrepare(filename = "../../data/mlml_output/taps_tapsbeta_mlml_output_filtered_1kbbin.tsv.gz")
taps[,"c"] <- 1-taps[,"hmc"]-taps[,"mc"]
taps <- na.omit(taps)
density <- TernaryDensity(taps, resolution = 25L, direction = 1)
pdf("plots/ternary_plot20200705/taps_tapsbeta_mlml.pdf")
TernaryPlot(axis.labels = seq(0, 100, by = 10), atip = "5hmC",btip = "C",ctip = "5mC")
ColourTernary(density, spectrum = col)
dev.off()

taps <- TernaryDiagramPrepare(filename = "../../data/mlml_output/tapsbeta_caps_mlml_output_filtered_1kbbin.tsv.gz")
taps[,"c"] <- 1-taps[,"hmc"]-taps[,"mc"]
density <- TernaryDensity(taps, resolution = 25L)
pdf("plots/ternary_plot20200705/tapsbeta_caps_mlml.pdf")
TernaryPlot(axis.labels = seq(0, 100, by = 20), atip = "5hmC",btip = "C",ctip = "5mC")
ColourTernary(density, spectrum = col)
dev.off()

taps <- TernaryDiagramPrepare(filename = "../../data/mlml_output/wgbs_ace_mlml_output_filtered_1kbbin.tsv.gz")
density <- TernaryDensity(taps, resolution = 25L)
pdf("plots/ternary_plot20200705/wgbs_ace_mlml.pdf")
TernaryPlot(axis.labels = seq(0, 100, by = 10), atip = "5hmC",btip = "C",ctip = "5mC")
ColourTernary(density, spectrum = col)
dev.off()

taps <- TernaryDiagramPrepare(filename = "../../data/mlml_output/wgbs_tabseq_mlml_output_filtered_1kbbin.tsv.gz")
density <- TernaryDensity(taps, resolution = 25L)
pdf("plots/ternary_plot20200705/wgbs_tabseq_mlml.pdf")
TernaryPlot(axis.labels = seq(0, 100, by = 10), atip = "5hmC",btip = "C",ctip = "5mC")
ColourTernary(density, spectrum = col)
dev.off()


taps <- TernaryDiagramPrepare(filename = "../../data/mlml_output/wgbs_oxbs_tabseq_mlml_output_filtered_1kbbin.tsv.gz")
density <- TernaryDensity(taps, resolution = 25L)
pdf("plots/ternary_plot20200705/wgbs_oxbs_tabseq_mlml.pdf")
TernaryPlot(axis.labels = seq(0, 100, by = 10), atip = "5hmC",btip = "C",ctip = "5mC")
ColourTernary(density, spectrum = col)
dev.off()


# pdf("plots/ternary_plot20200705/background.pdf")
# TernaryPlot(axis.labels = seq(0, 100, by = 10), atip = "5hmC",btip = "C",ctip = "5mC")
# dev.off()

################################################
################################################
######### below is obsolete ############


date  = format(Sys.time(), "%d%b%Y") # date
output  = "" # output file name
samples = c("TAPS", "YL711505m","YL712508hm") # sample names 
assays = c("TAPS", "TAPSbeta", "CAPS") # assays  
ncores = 4 # n CPUs
data_path = "../../data/" # data/
cutoff_reset <- data.frame(
  p = rep(c(1, 0.05, 0.01),3),
  depth = c(3,3,3, 5,5,5, 10,10,10), 
  mod_level = 0.1) # combinations of cutoffs
cutoff <- c(depth = 5)
fp_rate <- read.table("false_positive_table20200415.txt", header = T, stringsAsFactors = F)
fp_rate <- rbind(fp_rate, c("TAPS", 0.0023)) # use 0.23% as TAPS's false postive rate

## preprocessing
preprocess <- function(df){
  colnames(df) <- c("chr","start","end","mod_level","n_mod","n_unmod")
  
  df$mod_level <- as.numeric(df$mod_level) # mod_level (contained '*')
  df$n_reads <- df$n_mod + df$n_unmod
  df$n_reads[is.na(df$n_reads)] <- 0
  return(df)
}

## binomial testing
binom_t <- function(df, assay) {
  
  p <- as.numeric(fp_rate$conversion[fp_rate$methods == assay]) # false positive rate/probability
  df$p_value <- mcmapply(df$n_mod, df$n_reads, 
                         FUN = function(x, y) {
                           if(y > 0) {
                             return(binom.test(x = x, n = y, p = p, alternative = "greater")$p.value)
                           } else {
                             return(NA)
                           }
                         }, 
                         mc.cores = ncores)
  
  df$p_adjusted <- p.adjust(df$p_value, method = "BH")
  return(df)
}

## read in astair output
chr=2
df_mod <- list()
for(itor in 1:length(samples)) {
  fn <- paste(data_path, samples[itor], "_dedup_mCtoT_chr", chr, "_CpG_masked_columnSeleted.tsv", sep = "") # TAPS_beta file name
  df_mod[[itor]] <- fread(file = fn, nThread = ncores,  stringsAsFactors = F, data.table = F)
  df_mod[[itor]] <- preprocess(df_mod[[itor]]) # preprocess data
  if(local) {
    df_mod[[itor]] <- df_mod[[itor]][ df_mod[[itor]]$start < 1e7,] ## select 
  }
  df_mod[[itor]] <- binom_t(df = df_mod[[itor]], assay = assays[itor]) # binomial test
}
names(df_mod) <- assays
df <- do.call("rbind", lapply(df_mod, function(x) return(x[,c(1:3)])))
df <- unique(df)
df <- df[order(df$start),]

for(itor in 1:length(samples)) { 
  df_mod[[itor]] <- df_mod[[itor]][match(df$start,  df_mod[[itor]]$start),]
  df_mod[[itor]]$chr    <-  df$chr
  df_mod[[itor]]$start  <-  df$start
  df_mod[[itor]]$end    <-  df$end
  
  df_mod[[itor]]$mod_level[is.na(df_mod[[itor]]$mod_level)] <- 0
  df_mod[[itor]]$n_reads[is.na(df_mod[[itor]]$n_reads)] <- 0
  df_mod[[itor]]$p_adjusted[is.na(df_mod[[itor]]$p_adjusted)] <- 1
  
}

idx <-  df_mod[[1]]$n_reads >= 5 &
        df_mod[[2]]$n_reads >= 5 &
        df_mod[[3]]$n_reads >= 5 

for (itor in 1:length(samples)) {
  df_mod[[itor]] <- df_mod[[itor]][idx, ]
}

df <- df[idx,]

## TAPS - TAPSbeta
df_plot <- data.frame(c   = 1 - df_mod[[1]]$mod_level,
                      mc  = df_mod[[2]]$mod_level,
                      hmc = df_mod[[1]]$mod_level - df_mod[[2]]$mod_level)

df_plot$hmc[df_plot$hmc < 0] <- 0
df_plot <- as.data.frame(t(apply(df_plot, 1, function(x) return(x/sum(x)))))



df_plot <- subset(df_plot, hmc > 0)
df_plot <- subset(df_plot, mc > 0)
df_plot <- subset(df_plot, c > 0)
## Ternar
coordinates <- df_plot

library('Ternary')
par(mar = rep(0.2, 4))
TernaryPlot(axis.labels = seq(0, 10, by = 1), alab = 'C', blab = 'mC', clab = 'hmC')

# nPoints <- 4000L
# coordinates <- cbind(abs(rnorm(nPoints, 2, 3)),
#                      abs(rnorm(nPoints, 1, 1.5)),
#                      abs(rnorm(nPoints, 1, 0.5)))
# 

ColourTernary(TernaryDensity(coordinates, resolution = 30L))
# TernaryPoints(coordinates, col = 'red', pch = '.')
TernaryDensityContour(coordinates, resolution = 36L)
  
  
## ggtern
## https://www.geo.fu-berlin.de/en/v/soga/Introduction-to-R/Plotting-Data/ternary-diagrams/index.html






sessionInfo()
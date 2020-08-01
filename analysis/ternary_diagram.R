#!/usr/bin/env Rscript
## Tenary diagram 
## author: Zhiyuan Hu
## created 23 Apr 2020
## last modified 4 july 2020
## creat ternary plot


### ternary plot
## output: plots/ternary_plot20200705


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



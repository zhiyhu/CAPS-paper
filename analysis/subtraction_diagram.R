#!/usr/bin/env Rscript
## Tenary diagram 
## author: Zhiyuan Hu
## created 23 Apr 2020
## last modified 18 Nov 2020
## create scatterplot coloured by density

library(data.table)
library(ggplot2)
library(cowplot)

# read data function
readData <- function(filename) {
  taps <- fread(filename, 
                stringsAsFactors = F, data.table = F, header = F, sep = "\t", verbose = T, nThread = 4)
  colnames(taps) <- c('chr','start','end','mc','hmc','c')
  taps$c <- as.numeric(taps$c)
  taps <- taps[!is.na(taps$c),]
  taps$mc <- as.numeric(taps$mc)
  taps$hmc <- as.numeric(taps$hmc)
  # taps <- taps[taps$c >0 & taps$hmc >0  & taps$mc > 0, ]
  m <- as.matrix(taps[,c(5,6,4)]) # in order of hmc, c, mc
  m <- na.omit(m)
  m <- as.data.frame(m)
  m$density <- densCols(m$mc, m$hmc, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  return(m)
}

# legend
pdf("plots_revision/substraction20201108/color_gradient.pdf")
plot(rep(1,100),col=(colorRampPalette(rev(rainbow(10, end = 4/6)))(100)), pch=15,cex=2)
dev.off()

# scatterplot
taps <- readData(filename = "../../data/mlml_output/taps_tapsbeta_mlml_output_filtered_1kbbin.tsv.gz")
p1 <- ggplot(taps) + geom_point(aes(x = mc, y = hmc, col = density), size = 1) + scale_color_identity() +  theme_bw() + xlab("5mC raw signal") + ylab("5hmC raw signal") + labs(title = "TAPS and TAPSbeta (substraction)") + xlim(0, 1) + ylim(0, 0.5)
 
taps <- readData("../../data/mlml_output/tapsbeta_caps_mlml_output_filtered_1kbbin.tsv.gz")
p2 <- ggplot(taps) + geom_point(aes(x = mc, y = hmc, col = density), size = 1) + scale_color_identity() +  theme_bw() + xlab("5mC raw signal") + ylab("5hmC raw signal") + labs(title = "TAPS and CAPS")+ xlim(0, 1) + ylim(0, 0.5)

taps <- readData("../../data/mlml_output/wgbs_ace_mlml_output_filtered_1kbbin.tsv.gz")
p3 <- ggplot(taps) + geom_point(aes(x = mc, y = hmc, col = density), size = 1) + scale_color_identity() +  theme_bw() + xlab("5mC raw signal") + ylab("5hmC raw signal")+ labs(title = "WGBS and ACE-seq") + xlim(0, 1) + ylim(0, 0.5)

taps <- readData("../../data/mlml_output/wgbs_tabseq_mlml_output_filtered_1kbbin.tsv.gz")
p4 <- ggplot(taps) + geom_point(aes(x = mc, y = hmc, col = density), size = 1) + scale_color_identity() +  theme_bw() + xlab("5mC raw signal") + ylab("5hmC raw signal") + labs(title = "WGBS and TAB-seq") + xlim(0, 1) + ylim(0, 0.5)

plot_grid(p1,p2,p3,p4, nrow = 1)
ggsave("plots_revision/substraction20201108/mlml_substraction_density_scatterplot.pdf", width = 12, height = 3, useDingbats = FALSE)



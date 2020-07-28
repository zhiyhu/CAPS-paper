#!/usr/bin/env Rscript
## Plot the HPLC results (Figure S2)
## author: Zhiyuan Hu
## created 28 Jun 2020
## last modified 28 Jun 2020

library(ggplot2)
library(RColorBrewer)

my_cols <- c("#0072B2","#CC79A7","#D55E00","#009E73")

data  <- read.csv("plots/HPLC_20200628/HPLC_res.csv", as.is = T)
head(data )

data$Methylation <- factor(data$Methylation, levels = c("5mC","5hmC","5fC","5caC"))

data$Reaction <- factor(data$Reaction, levels = c("mESC_control","1x Oxidation","2x Oxidation"))

# display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE,
#                                        colorblindFriendly=F)
# mycol <- RColorBrewer::brewer.pal(12, "Paired")[c(2,12,10,4,7)]
# # mycol <- RColorBrewer::brewer.pal(8, "Set2")[c(1,2,3,4,6)]

ggplot(data = data, aes(x = Reaction, y = Percentage, fill = Methylation, col = Methylation)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  # geom_point(aes(fill = Methylation),  size = 1.2 , alpha = 0.7,
  #            shape = 21, position=position_jitterdodge()) +
  scale_fill_manual(values=my_cols)+scale_color_manual(values=my_cols)+
  geom_text(aes(label=round(Percentage, 1)), vjust=-0.4) +
  xlab("") +
  ylab("Species percentage (%)") + 
  theme_light() + theme(axis.text = element_text(color = "black"))

ggsave(filename="plots/HPLC_20200628/HPLC_barplot.pdf", useDingbats=FALSE, height=3,width=5.5)


df_plot <- data.frame(reaction = c("1x oxidation", "2x oxidation"),
                      conversion = c(82.81988, 97.16548))
ggplot(data = df_plot , aes(x = reaction, y = conversion)) + 
  geom_bar(stat = "identity", col = "black") +
  geom_text(aes(label=round(conversion, 1)), vjust=-0.4) +
  xlab("") +
  ylab("Conversion rate (%)") + 
  theme_classic() + theme(axis.text = element_text(color = "black"))

ggsave(filename="plots/HPLC_20200628/HPLC_barplot_conversion.pdf", useDingbats=FALSE, height=3.5,width=1.8)

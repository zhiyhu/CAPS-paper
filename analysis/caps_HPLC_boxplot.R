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


# 
# df2 <- read.csv("data/Figure 2 HPLC-MS_plot_SD.csv", as.is = T)
# head(df2)
# df2$Methylation <- factor(df2$Methylation, levels = c("5mC","5hmC","5fC","5caC","DHU"))
# 
# df2$Reaction <- gsub(df2$Reaction, pattern = "Tet", replacement = "NgTET1")
# df2$Reaction <- factor(df2$Reaction, levels = c("mESC_control","NgTET1","pic-borane"))
# 


# p <- ggplot(data = df2, aes(x = Reaction, y = Percentage, fill = Methylation)) + 
#     geom_bar(stat="identity", color="black", 
#              position=position_dodge()) +
#     geom_errorbar(aes(ymin=Percentage-SD, ymax=Percentage+SD), width=.4,
#                   position=position_dodge(.9), size=1) +
#     geom_point(data = data, aes(x = Reaction, y = Percentage, fill = Methylation),  size = 2 , alpha = 0.8,
#                shape = 21, position=position_jitterdodge()) +
#     scale_fill_manual(values=mycol)+
#     xlab("") +
#     ylab("Species percentage") + scale_y_continuous(labels = scales::percent) +
#     theme_classic() + 
#     theme(legend.text = element_text(size=15, face="bold", margin = margin(t = 10)),
#           # legend.key = element_rect(color = NA, fill = NA),
#           # legend.key.size = unit(0.8, "cm"),
#           legend.title = element_blank(),
#           legend.spacing.x = unit(0.2,"cm"),
#           # legend.spacing.y = unit(0.2, "cm"),
#           axis.text.x = element_text(size=15, face="bold", color = "black"),
#           axis.text.y = element_text(size=15, face="bold", color = "black"),
#           axis.title.y = element_text(size=15, face="bold"))
# 
# p
# ggsave(plot=p,height=4.5,width=8,dpi=300, filename="barplot.pdf", useDingbats=FALSE)

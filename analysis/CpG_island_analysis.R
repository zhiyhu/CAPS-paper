#!/usr/bin/env Rscript
## Analyse the modification levels and coverage at and around the CpG islands
## Zhiyuan Hu
## 20 May 2020


########### run CpG_island_prepare_bins.R (local) -----------------
########### run CpG_island_analysis.sh -----------------


############################  remote --------------------

library(data.table)
library(RColorBrewer)
library(tidyr)
library(ggplot2)

my_col <- brewer.pal(name = "Set2", n = 8)

sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
} else {
  os <- "No"
}

fn <- paste("../tidy_data/cpg_island/5hmC_caps_mlml_ace_tab_cpgIsland", c("_bins","_left_4kb_flank_20wins","_right_4kb_flank_20wins"), "_annotated_", sep = "")
chrs <- paste("chr",c(1:19, "X"), sep = "")


all <- c()

## chr     start   end     strand  caps.mod        caps.unmod      ace.mod ace.unmod       tab.mod tab.unmod       mlml.pm mlml.ph mlml.pu mlml.conflicts

# merge files from different chrs
for(itor in chrs){
  
  fn_tab  <- paste(fn, itor, ".bed.gz", sep = "") 
  
  tmp <- fread(fn_tab[1], header = F, nThread = 4,stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
  tmp$V18 <- paste("cpgisland_",  tmp$V18, sep = "")
  tab <- tmp
  
  tmp <- fread(fn_tab[2], header = F, nThread = 4,stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
  tmp$V18 <- paste("leftflank_",  tmp$V18, sep = "")
  tab <- rbind(tab, tmp[,c(1:19)])
  
  tmp <- fread(fn_tab[3], header = F, nThread = 4,stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
  tmp$V18 <- paste("rightflank_",  tmp$V18, sep = "")
  tab <- rbind(tab,  tmp[,c(1:19)])
  rm(tmp)
  
  all  <- rbind(all,  tab)
  
}

# calculate depth from mod and unmod reads
all$caps_depth <- all$V5 +all$V6
all$ace_depth <- all$V7 + all$V8
all$tab_depth <- all$V9 + all$V10

## df order for the summary table
df_order <- data.frame(region = c(rep("leftflank", 20), rep("cpgisland", 10), rep("rightflank", 20)),
                       bin = paste("bin", c(1:20, 1:10, 1:20), sep = ""),
                       number = 1:50)

df_order$name <- paste(df_order$region, df_order$bin, sep = "_")

all$caps_depth_downsampled <- all$caps_depth/sum(all$caps_depth) * sum(all$tab_depth)
all$ace_depth_downsampled <- all$ace_depth/sum(all$ace_depth) * sum(all$tab_depth)

## average coverage of caps
df_summary <- c()
tmp <- aggregate(caps_depth_downsampled~V18, all[!is.na(all$V18),], FUN = mean) # V16 is annotation; V12 is total depth
colnames(tmp) <- c("region","stat")
df_summary <- rbind(df_summary, tmp)
df_order$caps_cov <- df_summary$stat[match(df_order$name, df_summary$region)]

## average coverage of ace
df_summary <- c()
tmp <- aggregate(ace_depth_downsampled~V18, all[!is.na(all$V18),], FUN = mean)
colnames(tmp) <- c("region","stat")
df_summary <- rbind(df_summary, tmp)
df_order$ace_cov <- df_summary$stat[match(df_order$name, df_summary$region)]

## average coverage of tabseq
df_summary <- c()
tmp <- aggregate(tab_depth~V18, all[!is.na(all$V18),], FUN = mean)
colnames(tmp) <- c("region","stat")
df_summary <- rbind(df_summary, tmp)
df_order$tab_cov <- df_summary$stat[match(df_order$name, df_summary$region)]

## count per bin
tmp <- data.frame(table(all$V18))
df_order$n <- tmp$Freq[match(df_order$name, tmp$Var1)]


df_plot <- gather(data = df_order, assay, cov, caps_cov:tab_cov)
ggplot() + geom_vline(xintercept = 50.5, col = "grey")+ geom_vline(xintercept = 70.5, col = "grey")+ 
  geom_line(aes(x = number, y = cov, col = assay),
            data = df_plot, stat="identity")  + 
  scale_color_manual(values = my_col[c(5,2,3)]) + xlim(0, 50)+ 
  theme_linedraw() + 
  theme(text = element_text(size = 12), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "top")  

if(os != "Darwin") {
  ggsave("plots/cpg_island20200704/cpg_island_coverage_caps_vs_ace_tabseq.pdf", width = 8, height = 4)
  write.table(df_order,"plots/cpg_island20200704/cpg_island_df_order_plotdata.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
}



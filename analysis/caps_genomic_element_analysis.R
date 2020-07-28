#!/usr/bin/env Rscript
## colocolization analysis with genomic regulatory elements of 5hmC
## Zhiyuan Hu
## 20 May 2020
## last modified: 26 July 2020
## ran remotely

library(data.table)
library(parallel)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
library(RColorBrewer)

depth_cutoff <- 5
p_cutoff <- 0.05
tabseq_fp <- 0.0222 ## TAB-seq false positive rate
ace_fp    <- 0.0047 ## ACE-seq false positive rate
caps_fp   <- 0.00716 # CAPS CpG false positive rate
ncores <- 4

## set up ----
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
} else {
  os <- "No"
}

## part 1 -----
## preprocess HMM file before R/genomic_element_analysis.sh
# df_element <- fread("../../data/mESC_cStates_HMM", skip = 1, sep = "\t", nThread = 4, data.table = F)
# fwrite(df_element <- df_element[,1:4], "../../data/mESC_cStates_HMM.bed.gz", col.names = F, row.names = F, sep = "\t")
# /home/obgynae/zyhu/projects/taps/ref/mESC_cStates_HMM_sorted.bed.gz

## run R/genomic_element_analysis.sh ---- 

## part 2 ----

if(os == "Darwin") {
  fn_caps <- "../../data/methy_data/caps_CpG_mods_annotated_chr2_top100k.bed.gz"
} else {
  fn_caps <- "../methy_data/caps_CpG_mods_annotated.bed.gz"
}

caps_all <- fread(fn_caps, header = F, nThread = 4, stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
colnames(caps_all) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth","chr_hmm","start_hmm","end_hmm","status","no")

## filter variants and depth
caps_all <- caps_all[caps_all$snv == "No",] ## filter by SNV
caps_all$depth <- caps_all$mod + caps_all$unmod
caps_all <- caps_all[caps_all$depth >= depth_cutoff,] # filter by depth
caps_all <- caps_all[,c("chr","start","end","mod_level","mod","unmod","depth","start_hmm","end_hmm","status","no")]

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

print("binomial test caps")
caps_all$p <- binom_t(depth = caps_all$depth, mod = caps_all$mod, p = caps_fp)
caps_all$fdr <- p.adjust(caps_all$p, method = "BH")
caps_all$status <- gsub('[0-9]+_', '', caps_all$status)# remove the number IDs from status names

fwrite(caps_all, "../tidy_data/5hmC_caps_hmm_annotated_binomial.bed.gz", row.names = F, col.names = F, sep = "\t")

## ## ## ## ## ## ## ## ## #### ## ## ## ## ## ## ## ## ##
## Count by genomic elements ------
## ## ## ## ## ## ## ## ## #### ## ## ## ## ## ## ## ## ##

## read caps data
fn_caps <- "../tidy_data/5hmC_caps_hmm_annotated_binomial.bed.gz"
caps_all <- fread(fn_caps, header = F, nThread = 4, stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
colnames(caps_all) <- c("chr","start","end","mod_level","mod","unmod","depth","start_hmm","end_hmm","status","no","p","fdr")
# caps_all <- caps_all[!is.na(caps_all$chr),]

head(caps_all)

##########  Count 5hmC sites -----

## caluclate number of 5hmCs observed
n <- sum(caps_all$fdr <= p_cutoff)

print(n)

## calculate the seq depth/coverage in each element type
df_tmp <- caps_all[caps_all$fdr <= p_cutoff,]

# df_element <- aggregate(depth~status, data = caps_all, FUN = sum) # coverage
hmm <- caps_all[,c("chr","start_hmm","end_hmm","status")]
hmm <- distinct(hmm)
hmm$length <- hmm$end_hmm - hmm$start_hmm

df_element <- aggregate(length~status, data = hmm, FUN = sum)

tmp <- as.data.frame(table(df_tmp$status))
df_element$n_sites <- tmp$Freq[match(df_element$status, tmp$Var1)] #
df_element$sites_perkb <- df_element$n_sites/df_element$length * 1e6

## length for each element in the whole genomes

# tmp <- table(caps_all$status)
# df_element$caps_cpg_covered <- tmp[match(df_element$status, names(tmp))]
# tmp <- table(caps_all$status[caps_all$caps.p <= p_cutoff])
# df_element$caps_cpg_5hmc <- tmp[match(df_element$status, names(tmp))]

##########  10 random sampling ----
set.seed(123545)
df_sampling <- replicate(10,sample(1:nrow(caps_all), size = n, replace = F))

sampling_df <- c()

for(itor in 1:10){
  tmp <- table(caps_all$status[df_sampling[,itor]])
  tmp <- as.numeric(tmp[match(df_element$status, names(tmp))])
  tmp <- tmp/df_element$length * 1e6
  sampling_df <- cbind(sampling_df, tmp)
}

df_element$sampling_freq <- rowMeans(sampling_df, na.rm = T)
df_element$sampling_err <- apply(sampling_df, 1, sd)

df_element$ratio <- df_element$sites_perkb/df_element$sampling_freq

if(os != "Darwin") {
  write.table(df_element, "plots/genomic_element20200628/df_element_summary.txt", row.names = F, col.names = T, sep = "\t", quote = F )
}

############# bar plot ----

# df_element <- read.delim("plots/genomic_element20200628/df_element_summary.txt", as.is = T)

df_element$sites_perkb_sd <- 0

df_plot <- data.frame(status = rep(df_element$status, 2),
                      n  = c(df_element$sites_perkb, df_element$sampling_freq),
                      sd = c(df_element$sites_perkb_sd, df_element$sampling_err),
                      group = rep(c("Observed", "Random (10x)"), each = 10),
                      ratio = c(rep("", 10), round(df_element$ratio, 2)))

df_plot$status <- factor(df_plot$status, levels = c("Active_Promoter", "Poised_Promoter","Strong_Enhancer","Poised_Enhancer","Txn_Transition","Txn_Elongation","Weak_Txn","Insulator","Repressed","Heterochrom"))
## bar plot
my_col <- brewer.pal(name = "Set2", n = 8)

ggplot(df_plot, aes(x = status, y = n, fill = group)) + 
  geom_bar(stat="identity", position=position_dodge(), col = "black") + 
  geom_errorbar(aes(ymin=n-sd, ymax=n+sd), width=.2,
                position=position_dodge(.9)) +
  theme_light() + 
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8, size = 12,colour = "black"),
        axis.text.y = element_text(colour = "black"))  + 
  geom_text(aes(label = ratio, x = status, y = max(n)), vjust = -0.5) + 
  scale_fill_manual(values = c("grey40","white")) + 
  xlab("") + ylab("Number of called high-confidence 5hmCGs\nin genomic feature")

ggsave("plots/genomic_element20200628/genomic_element_5hmC_vs_random_barplot.pdf", width = 10, height = 6)

if(os != "Darwin") {
  write.table(df_plot, "plots/genomic_element20200628/genomic_element_5hmC_vs_random_barplot_data.txt", row.names = F, col.names = T, sep = "\t", quote = F )
}

######## pie plot -------

df_summary <- table(caps_all$status[caps_all$fdr <= p_cutoff])
df_summary <- data.frame(df_summary)

df_summary$Var1 <- factor(df_summary$Var1, levels = rev(c("Active_Promoter", "Poised_Promoter","Strong_Enhancer","Poised_Enhancer","Txn_Transition","Txn_Elongation","Weak_Txn","Insulator","Repressed","Heterochrom")))
# df_summary <- df_summary[!df_summary$Var1 %in% c("Heterochrom"),]

df_summary$percent <- df_summary$Freq/sum(df_summary$Freq)
df_summary <- df_summary[match(c("Active_Promoter", "Poised_Promoter","Strong_Enhancer","Poised_Enhancer","Txn_Transition","Txn_Elongation","Weak_Txn","Insulator","Repressed","Heterochrom"), df_summary$Var1),]

ggplot(df_summary, aes(x="", y=percent, fill=Var1))+ coord_polar("y", start=0)+
  geom_bar(width = 1, stat = "identity") + scale_fill_brewer(palette="Paired") + theme_minimal() + xlab("") + ylab("") +
  theme(text = element_text(size = 12, colour = "black"))

ggsave("plots/genomic_element20200628/genomic_element_5hmC_piechart.pdf", width = 10, height = 6)


if(os != "Darwin") {
  write.table(df_summary, "plots/genomic_element20200628/genomic_element_5hmC_piechart_data.txt", row.names = F, col.names = T, sep = "\t", quote = F )
}



#### boxplot: distribution

# caps_all$mod_level <- as.numeric(caps_all$mod_level)
# df_quantile <- aggregate(mod_level~status, data = caps_all[caps_all$caps.p <= p_cutoff,], FUN = summary)
# df_quantile <- as.data.frame(df_quantile)
# df_summary <- df_quantile$mod_level 
# colnames(df_summary) <- gsub(pattern = "[.]", replacement = "", x = colnames(df_summary))
# df_summary <- data.frame(df_summary)
# df_summary$status <- df_quantile$status
# df_summary$status <- factor(df_summary$status, levels = c("Active_Promoter", "Poised_Promoter","Strong_Enhancer","Poised_Enhancer","Txn_Transition","Txn_Elongation","Weak_Txn","Insulator","Repressed","Heterochrom"))
# 
# ggplot(df_summary, aes(x = status, ymin = Min, lower = X1st.Qu, middle = Median, upper = X3rd.Qu, ymax = Max)) +
#   geom_boxplot(stat = "identity") + theme_linedraw() +
#   theme(text = element_text(size = 12, colour = "black"),
#         axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8, size = 12))   + 
#   xlab("") + ylab("Modification level of high-confidence CpGs")
# ggsave("plots/genomic_element20200520/genomic_element_5hmC_distribution_boxplot.pdf", width = 10, height = 6)
# 
# 
# if(os != "Darwin") {
#   write.table(df_summary, "plots/genomic_element20200520/genomic_element_5hmC_distribution_boxplot_data.txt", row.names = F, col.names = T, sep = "\t", quote = F )
# }


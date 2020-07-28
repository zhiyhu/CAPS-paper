#!/usr/bin/env Rscript
## preprocessing CAPS data, ACE data and TABseq data
## Zhiyuan Hu
## 11 May 2020
## last modified 14 May 2020

library(data.table)
library(parallel)

depth_cutoff <- 0
p_cutoff <- 0.001
tabseq_fp <- 0.0222 ## TAB-seq false positive rate
ace_fp <- 0.0047 ## ACE-seq false positive rate
ncores <- 8

## set up ----
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
} else {
  os <- "No"
}

if(os == "Darwin") { ## if run locally (mock run only with chr 2)
  fn_cpg  <- "../../data/mm9_newspike_CpG_chr2.bed.gz"
  fn_caps <- "../../data/caps/YL712508hm_dedup_chrselected_mCtoT_CpG_masked_chr2.mods.gz"
  fn_ace  <- "../../data/ace/GSE116016_ACE-Seq_WT.ESC.mm9_CG_chr2.txt.gz"
  fn_mlml <- "../../data/mlml_output/taps_tapsbeta_mlml_output_chr2.tsv.gz"
  fn_taps <- "../../data/data_for_mlml/tidy_taps_mlml_chr2.tsv.gz"
  fn_beta <- "../../data/data_for_mlml/tidy_tapsbeta_mlml_chr2.tsv.gz"
  fn_tab  <- "../../data/tabseq/tabseq_bismark_bt2_dedup2_sorted.bismark_chr2.cov.gz"
  chrs <- "chr2"
} else {
  fn_cpg  <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/mm9_CpG_sites.bed.gz"
  fn_caps <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/caps_CpG_mods.bed.gz"
  fn_ace  <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/ace_CpG_mods.bed.gz"
  fn_mlml <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/taps_tapsbeta_mlml_CpG_mods.bed.gz"
  fn_taps <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/taps_CpG_mods.bed.gz"
  fn_beta <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/tapsbeta_CpG_mods.bed.gz"
  fn_tab  <- "/home/obgynae/zyhu/projects/taps/analysis_caps/methy_data/tabseq_CpG_mods.bed.gz"
chrs <- paste("chr", c(1:19, "X","Y"), sep = "")
}


##################### part 1 -----------------

## CpG list ----
print("reading cpg")
cpg <- fread(fn_cpg, stringsAsFactors = F, data.table = F, header = F, sep = "\t", nThread = 4, verbose = T)
colnames(cpg) <- c("chr", "start", "end", "strand","full_context","context")

cpg_all <- cpg
rm(cpg)

## CAPS data ----
print("reading caps")

caps_all <- fread(fn_caps, stringsAsFactors = F, data.table = F, header = F, sep = "\t", nThread = 4, verbose = T)
colnames(caps_all) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")

tmp <- c()

for(itor in chrs){
  print(itor)
  cpg  <- cpg_all [cpg_all $chr == itor, ]
  caps <- caps_all[caps_all$chr == itor, ]
  
  ## caps
  cpg$caps.mod <- caps$mod[match(cpg$start, caps$start)]
  cpg$caps.unmod <- caps$unmod[match(cpg$start, caps$start)]
  cpg$snv <- caps$snv[match(cpg$start, caps$start)] ## add SNV info
  
  sum(caps$start %in% cpg$start)
  
  ## remove SNV
  cpg <- cpg[!is.na(cpg$caps.mod),]
  cpg <- cpg[cpg$snv == "No",]
  cpg <- cpg[,c(1:4,7,8)]
  tmp <- rbind(tmp, cpg)
}

cpg_all <- tmp
# cpg_all$caps.depth <- cpg_all$caps.mod + cpg_all$caps.unmod

rm(caps, caps_all, tmp)
gc()

## ACE data ----
print("reading ACE data")

ace_all <- fread(fn_ace, stringsAsFactors = F, data.table = F, header = T, sep = "\t", nThread = 4, verbose = T)
colnames(ace_all) <- c("chr","start","end","mod_level","depth","strand")
ace_all$depth <- as.numeric(ace_all$depth)
ace_all$mod_level <- as.numeric(ace_all$mod_level)
ace_all$mod   <- round(ace_all$depth * ace_all$mod_level, digits = 0)
ace_all$unmod <- ace_all$depth - ace_all$mod

tmp <- c()

for(itor in chrs){
  print(itor)
  cpg <- cpg_all[cpg_all$chr == itor, ]
  ace <- ace_all[ace_all$chr == itor, ]
  
  ## ace
  cpg$ace.mod <- ace$mod[match(cpg$start, ace$start)]
  cpg$ace.unmod <- ace$unmod[match(cpg$start, ace$start)]

  sum(ace$start %in% cpg$start)
  tmp <- rbind(tmp, cpg)
}

cpg_all <- tmp
# cpg_all$ace.depth  <- cpg_all$ace.mod  + cpg_all$ace.unmod

rm(ace_all, ace, tmp)
gc()


## TAB-seq data ----
print("reading ACE data")

ace_all <- fread(fn_tab, stringsAsFactors = F, data.table = F, header = F, sep = "\t", nThread = 4, verbose = T)
colnames(ace_all) <- c("chr","start","end","mod_level","mod","unmod")
ace_all$start <- ace_all$start - 1

tmp <- c()

for(itor in chrs){
  print(itor)
  cpg <- cpg_all[cpg_all$chr == itor, ]
  ace <- ace_all[ace_all$chr == itor, ]
  
  ## ace
  cpg$tab.mod <- ace$mod[match(cpg$start, ace$start)]
  cpg$tab.unmod <- ace$unmod[match(cpg$start, ace$start)]
  
  sum(ace$start %in% cpg$start)
  tmp <- rbind(tmp, cpg)
}

cpg_all <- tmp
# cpg_all$tab.depth  <- cpg_all$tab.mod  + cpg_all$tab.unmod

rm(ace_all, ace, tmp)
gc()

## mlml data ----
print("reading mlml data")

ace_all <- fread(fn_mlml, stringsAsFactors = F, data.table = F, header = T, sep = "\t", nThread = 4, verbose = T)
colnames(ace_all) <- c("chr","start","end","pm","ph","pu","conflicts")

tmp <- c()

for(itor in chrs){
  print(itor)
  cpg <- cpg_all[cpg_all$chr == itor, ]
  ace <- ace_all[ace_all$chr == itor, ]
  
  ## ace
  cpg$mlml.pm <- ace$pm[match(cpg$start, ace$start)]
  cpg$mlml.ph <- ace$ph[match(cpg$start, ace$start)]
  cpg$mlml.pu <- ace$pu[match(cpg$start, ace$start)]
  cpg$mlml.conflicts <- ace$conflicts[match(cpg$start, ace$start)]

  sum(ace$start %in% cpg$start)
  tmp <- rbind(tmp, cpg)
}

cpg_all <- tmp
# cpg_all$tab.depth  <- cpg_all$tab.mod  + cpg_all$tab.unmod

rm(ace_all, ace, tmp)


head(cpg_all)


for(i in c(5:13)){
  cpg_all[is.na(cpg_all[,i]),i] <- 0
}

cpg_all$mlml.conflicts[is.na(cpg_all$mlml.conflicts)] <- 5

head(cpg_all)


if(os != "Darwin") {
  fwrite(cpg_all, "../tidy_data/5hmC_caps_mlml_ace_tab_mods.unfiltered.bed.gz", sep = "\t", row.names = F, col.names = T)
}



## binomial test (disabled) -----
## don't need 

# fp_rate <- read.table("false_positive_table.txt", header = T, stringsAsFactors = F)
# fp_rate <- rbind(fp_rate, c("TAPS", 0.0023)) # use 0.23% as TAPS's false postive rate
# 
# binom_t <- function(depth, mod, p) { # false positive rate/probability
#   p_value <- mcmapply(as.numeric(mod), as.numeric(depth), 
#                          FUN = function(x, y) {
#                            if(!is.na(x) & !is.na(y)) {
#                              if(y > 0) {
#                                return(binom.test(x = x, n = y, p = p, alternative = "greater")$p.value)
#                              } else {
#                                return(1)
#                              }
#                            } else {
#                              return(1)
#                            }
#                            
#                          }, 
#                          mc.cores = ncores)
#   
#   # p_adjusted <- p.adjust(p_value, method = "BH")
#   return(p_value)
# }

print("sum(is.na(cpg_all$chr)) = ")
sum(is.na(cpg_all$chr))

cpg_all <- cpg_all[!is.na(cpg_all$chr),] ## remove NA rows

# print("binomial test caps")
# caps_fp <- as.numeric(fp_rate$conversion[fp_rate$methods == "CAPS"])
# cpg_all$caps.p <- binom_t(depth = (cpg_all$caps.mod + cpg_all$caps.unmod), mod = cpg_all$caps.mod, p = caps_fp)
# print("binomial test ace")
# cpg_all$ace.p  <- binom_t(depth = (cpg_all$ace.mod + cpg_all$ace.unmod),  mod = cpg_all$ace.mod,  p = ace_fp)
# print("binomial test tabseq")
# cpg_all$tab.p  <- binom_t(depth = (cpg_all$tab.mod + cpg_all$tab.unmod),  mod = cpg_all$tab.mod,  p = tabseq_fp)

cpg_all$caps.p <- 0
cpg_all$ace.p  <- 0
cpg_all$tab.p  <- 0

if(os != "Darwin") {
  fwrite(cpg_all, "../tidy_data/5hmC_caps_mlml_ace_tab_mods_binomial_tested.unfiltered.bed.gz", sep = "\t", row.names = F, col.names = T)
} else {
  fwrite(cpg_all, "../../data/tidy_data/5hmC_caps_mlml_ace_tab_mods_binomial_tested.unfiltered_chr2.bed.gz", sep = "\t", row.names = F, col.names = F)
}


####################### part 2: clean up the columns ----

cpg_all <- fread("../tidy_data/5hmC_caps_mlml_ace_tab_mods_binomial_tested.unfiltered.bed.gz", header = T, nThread = 4,stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)

tidy_data <- cpg_all

tidy_data$caps.mod2 <- tidy_data$caps.mod
tidy_data$ace.mod2 <- tidy_data$ace.mod
tidy_data$tab.mod2 <- tidy_data$tab.mod


## remove mlml results with conflicts over 0
tidy_data$mlml.mod <- tidy_data$mlml.ph
tidy_data$mlml.mod[tidy_data$mlml.conflicts >= 1] <- 0
tidy_data$mlml.unmod <- 1 - tidy_data$mlml.ph
tidy_data$mlml.unmod[tidy_data$mlml.conflicts >= 1] <- 0


## remove binomial with p > 0.01
tidy_data$caps.mod[tidy_data$caps.p > p_cutoff] <- 0
tidy_data$ace.mod[tidy_data$ace.p > p_cutoff] <- 0
tidy_data$tab.mod[tidy_data$tab.p > p_cutoff] <- 0

if(os != "Darwin") {
  fwrite(tidy_data, "../tidy_data/5hmC_caps_mlml_ace_tab_mods.bed.gz", sep = "\t", row.names = F, col.names = T)
}
















################################################################################################
##################################### below is obsolete -------
################################################################################################


## caps vs tabseq vs ace -----

idx <- which(cpg_all$caps.depth >= depth_cutoff &
               cpg_all$ace.depth >= depth_cutoff &
               cpg_all$tab.depth >= depth_cutoff)

tidy_data <- cpg_all[idx, c("chr","start","end", "strand", "caps.mod", "caps.unmod", 
                            "ace.mod", "ace.unmod","tab.mod" ,"tab.unmod", "caps.p", "ace.p","tab.p")]

tidy_data$caps.mod2 <- tidy_data$caps.mod
tidy_data$ace.mod2 <- tidy_data$ace.mod
tidy_data$tab.mod2 <- tidy_data$tab.mod

tidy_data$caps.mod[tidy_data$caps.p > 0.05] <- 0
tidy_data$ace.mod[tidy_data$ace.p > 0.05] <- 0
tidy_data$tab.mod[tidy_data$tab.p > 0.05] <- 0

if(os != "Darwin") {
  fwrite(tidy_data, "../tidy_data/5hmC_caps_ace_tab_mods.bed.gz", sep = "\t", row.names = F, col.names = F)
}



## mlml vs tabseq vs ace -----


idx <- which(#cpg_all$mlml.depth >= depth_cutoff &
               cpg_all$caps.depth >= depth_cutoff&
               cpg_all$ace.depth >= depth_cutoff &
               cpg_all$tab.depth >= depth_cutoff & 
               cpg_all$mlml.conflicts == 0 & 
               !is.na(cpg_all$mlml.conflicts) &
               cpg_all$mlml.ph != "" & 
               !is.na(cpg_all$ mlml.ph))

tidy_data <- cpg_all[idx, c("chr","start","end", "strand", "mlml.mod", "mlml.unmod", 
                            "ace.mod", "ace.unmod","tab.mod" ,"tab.unmod", "ace.p","tab.p")]

tidy_data$ace.mod2 <- tidy_data$ace.mod
tidy_data$tab.mod2 <- tidy_data$tab.mod

tidy_data$ace.mod[tidy_data$ace.p > 0.05] <- 0
tidy_data$tab.mod[tidy_data$tab.p > 0.05] <- 0

if(os != "Darwin") {
  fwrite(tidy_data, "../tidy_data/5hmC_tapsmlml_ace_tab_mods.bed.gz", sep = "\t", row.names = F, col.names = F)
}

### all four

# depth_cutoff <- 3

idx <- which(#cpg_all$mlml.depth >= depth_cutoff &
    cpg_all$ace.depth >= depth_cutoff &
    cpg_all$tab.depth >= depth_cutoff & 
    cpg_all$mlml.conflicts == 0 & 
    !is.na(cpg_all$mlml.conflicts) &
    cpg_all$mlml.mod != "" & 
    !is.na(cpg_all$mlml.mod))


cpg_all$mlml.mod[is.na(cpg_all$mlml.mod)]

tidy_data <- cpg_all[idx, c("chr","start","end", "strand", "caps.mod", "caps.unmod","mlml.mod", "mlml.unmod", 
                            "ace.mod", "ace.unmod","tab.mod" ,"tab.unmod","caps.p", "ace.p","tab.p")]

tidy_data$caps.mod2 <- tidy_data$caps.mod
tidy_data$ace.mod2 <- tidy_data$ace.mod
tidy_data$tab.mod2 <- tidy_data$tab.mod

tidy_data$mlml.mod <- tidy_data$mlml.ph 
tidy_data$mlml.unmod <- 1 - tidy_data$mlml.ph
tidy_data$caps.mod[tidy_data$caps.p > 0.05] <- 0
tidy_data$ace.mod[tidy_data$ace.p > 0.05] <- 0
tidy_data$tab.mod[tidy_data$tab.p > 0.05] <- 0

if(os != "Darwin") {
  fwrite(tidy_data, "../tidy_data/5hmC_caps_tapsmlml_ace_tab_mods.bed.gz", sep = "\t", row.names = F, col.names = F)
}


sessionInfo()


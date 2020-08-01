#!/usr/bin/env Rscript
## preprocess and benchmark TAPSbeta and oxBS
## Zhiyuan Hu
## 1 May 2020 (last modified 27 June 2020)
## Also see benchmark_CAPS_TABseq.R
## run this script by run_5mc_benchmark.sh

library(data.table)

chrs <- paste("chr", c(1:19, "X","Y"), sep = "") # chrs
dir <- "$HOME/projects/taps/"
  
fn_cpg <-  paste(dir, "analysis_caps/methy_data/mm9_CpG_sites.bed.gz", sep = "") ## cpg sites
fn_beta <- paste(dir, "analysis_caps/methy_data/tapsbeta_CpG_mods.bed.gz", sep = "") # TAPS-beta file
fn_oxbs <- paste(dir, "analysis_caps/methy_data/oxbs_",
                   c(1:3), "_CpG_mods.bed.gz", sep = "")

fn <- c(fn_cpg)
fn_out <- "$HOME/projects/taps/analysis_caps/tidy_data/5mC_tapsbeta_oxbs_mods.bed.gz"



## CpG list ----
print("read cpg")
cpg <- fread(fn[1], stringsAsFactors = F, data.table = F, header = F, sep = "\t")
colnames(cpg) <- c("chr", "start", "end", "strand","full_context","context")

cpg_combined <- data.frame(chr = cpg$chr,
                           start = cpg$start,
                           end = cpg$end + 1)

cpg_combined <- cpg_combined[cpg_combined$end %in% cpg$end, ]
cpg_combined$strand1 <- cpg$strand[match(cpg_combined$start, cpg$start)] # + is C, - is G
cpg_combined$strand2 <- cpg$strand[match(cpg_combined$end, cpg$end)]
cpg_combined <- cpg_combined[paste(cpg_combined$strand1, cpg_combined$strand2) == "+ -",]

cpg_combined$start1 <- cpg_combined$start
cpg_combined$start2 <- cpg_combined$start + 1


## TAPS-beta data ----
print("read tapsbeta")
tapsb <- fread(fn_beta, stringsAsFactors = F, data.table = F, header = F, sep = "\t")
colnames(tapsb) <- c("chr","start","end","mod_level","mod","unmod","ref","mut","context","CpG","snv","depth")


split_column <- function(x, n) {  # n: the column
  if(is.na(x)) {
    return(0)
  } else {
    tmp <- unlist(strsplit(x = x, split = ","))[n]
    return(tmp)
  }
}


## oxBS data ----
print("read oxbs")
oxbs1 <- fread(fn_oxbs[1], stringsAsFactors = F, data.table = F, header = F) ## read in three files
oxbs2 <- fread(fn_oxbs[2], stringsAsFactors = F, data.table = F, header = F)
oxbs3 <- fread(fn_oxbs[3], stringsAsFactors = F, data.table = F, header = F)
colnames(oxbs1) <- colnames(oxbs2) <- colnames(oxbs3) <- c("chr","start","end","mod_level","mod","unmod")

## merge three outputs to achieve higher depth
oxbs <- merge(oxbs1, oxbs2, by = c("chr","start","end"), all = T)
oxbs <- merge(oxbs, oxbs3, by = c("chr","start","end"), all = T)

rm(oxbs1, oxbs2, oxbs3)

oxbs$mod.x[is.na(oxbs$mod.x)]     <- 0
oxbs$unmod.x[is.na(oxbs$unmod.x)] <- 0
oxbs$mod.y[is.na(oxbs$mod.y)]     <- 0
oxbs$unmod.y[is.na(oxbs$unmod.y)] <- 0
oxbs$mod[is.na(oxbs$mod)]         <- 0
oxbs$unmod[is.na(oxbs$unmod)]     <- 0

oxbs$mod_sum   <- oxbs$mod.x + oxbs$mod.y + oxbs$mod
oxbs$unmod_sum <- oxbs$unmod.x + oxbs$unmod.y + oxbs$unmod
oxbs <- oxbs[,c(1:3, 13:14)]
oxbs$start <- oxbs$start - 1 # adjust the range to match tapsbeta coordinates

sum(oxbs$start %in% cpg$start)
# [1] 1560612

## merge data ----
print("merge data")
cpg_all    <- cpg
tapsb_all  <- tapsb   
oxbs_all   <- oxbs   
cpg_new    <- c()

for(itor in chrs){
  print(itor)
  cpg     <- cpg_all    [cpg_all   $chr == itor, ]
  tapsb   <- tapsb_all  [tapsb_all $chr == itor, ]
  oxbs    <- oxbs_all   [oxbs_all  $chr == itor, ]
  
  ## taps-beta
  cpg$tapsb.mod <- tapsb$mod[match(cpg$start, tapsb$start)]
  cpg$tapsb.unmod <- tapsb$unmod[match(cpg$start, tapsb$start)]
  cpg$snv <- tapsb$snv[match(cpg$start, tapsb$start)] ## add SNV info
  
  sum(tapsb$start %in% cpg$start)
  # [1] 2899799
  
  ## remove SNV
  # cpg <- cpg[cpg$snv == "No",]
  
  ## oxbs
  cpg$oxbs.mod   <- oxbs$mod[  match(cpg$start, oxbs$start)]
  cpg$oxbs.unmod <- oxbs$unmod[match(cpg$start, oxbs$start)]
  
  cpg$tapsb.mod  [is.na(cpg$tapsb.mod)] <- 0
  cpg$tapsb.unmod[is.na(cpg$tapsb.unmod)] <- 0
  
  cpg$oxbs.mod  [is.na(cpg$oxbs.mod)] <- 0
  cpg$oxbs.unmod[is.na(cpg$oxbs.unmod)] <- 0
  
  cpg_new <- rbind(cpg_new, cpg)
  
}

rm(oxbs_all, oxbs, tapsb, tapsb_all, cpg_all, cpg)

print("merge C and G in CpGs")

merged_new <- c()

for(itor in chrs){

  print(itor)
  merged <- cpg_combined[cpg_combined$chr == itor, ]
  cpg <- cpg_new[cpg_new$chr == itor,]
  
  ## beta mod
  merged$beta.mod1 <- cpg$tapsb.mod[match(merged$start1, cpg$start)]
  merged$beta.mod2 <- cpg$tapsb.mod[match(merged$start2, cpg$start)]
  
  merged$beta.unmod1 <- cpg$tapsb.unmod[match(merged$start1, cpg$start)]
  merged$beta.unmod2 <- cpg$tapsb.unmod[match(merged$start2, cpg$start)]
  
  ## oxbs mod
  merged$oxbs.mod1 <- cpg$oxbs.mod[match(merged$start1, cpg$start)]
  merged$oxbs.mod2 <- cpg$oxbs.mod[match(merged$start2, cpg$start)]
  
  merged$oxbs.unmod1 <- cpg$oxbs.unmod[match(merged$start1, cpg$start)]
  merged$oxbs.unmod2 <- cpg$oxbs.unmod[match(merged$start2, cpg$start)]
  
  
  ## snv
  merged$snv1 <- cpg$snv[match(merged$start1, cpg$start)]
  merged$snv2 <- cpg$snv[match(merged$start2, cpg$start)]
  
  merged$beta.mod <- merged$beta.mod1 + merged$beta.mod2
  merged$beta.unmod <- merged$beta.unmod1 + merged$beta.unmod2
  merged$oxbs.mod <-   merged$oxbs.mod1   + merged$oxbs.mod2
  merged$oxbs.unmod <- merged$oxbs.unmod1 + merged$oxbs.unmod2
  
  merged$beta.mod[is.na(merged$beta.mod)] <- 0
  merged$oxbs.mod[is.na(merged$oxbs.mod)] <- 0
  merged$beta.unmod[is.na(merged$beta.unmod)] <- 0
  merged$oxbs.unmod[is.na(merged$oxbs.unmod)] <- 0
  
  ## remove mutated sites
  merged <- merged[which(merged$snv1 == "No" & merged$snv2 == "No"),]
  merged <- merged[,c("chr","start","end","beta.mod","beta.unmod","oxbs.mod","oxbs.unmod")]
  merged_new <- rbind(merged_new, merged)
}

fwrite(cpg_new, "../tidy_data/5mC_tapsbeta_oxbs_unfiltered_CpG_mods.bed.gz", col.names = T, sep = "\t")
fwrite(merged_new, "../tidy_data/5mC_tapsbeta_oxbs_unfiltered_CpGmerged_mods.bed.gz", col.names = T, sep = "\t")

################ benchmarking analysis #########################

rm(cpg, cpg_combined, cpg_new, merged)
gc()

## depth = 10
depth_cutoff <- 10

merged_keep <- merged_new

merged_new <- merged_new[(merged_new$beta.mod + merged_new$beta.unmod) >= depth_cutoff & 
                           (merged_new$oxbs.mod + merged_new$oxbs.unmod) >= depth_cutoff,]

merged_new$beta_level <- merged_new$beta.mod/(merged_new$beta.mod + merged_new$beta.unmod)
merged_new$oxbs_level <- merged_new$oxbs.mod/(merged_new$oxbs.mod + merged_new$oxbs.unmod)

df_res <- data.frame(assay1 = "oxBS",
                     assay2 = "TAPS-beta",
                     pearson_cor = cor(merged_new$beta_level, merged_new$oxbs_level),
                     spearman_cor = cor(merged_new$beta_level, merged_new$oxbs_level, method = "spearman"),
                     CpG_sites_used = nrow(merged_new))

write.csv(df_res, "plots/5mc_tapsbeta_benchmark20200627/correlation_result_depth10.csv")

pdf("plots/5mc_tapsbeta_benchmark20200627/smoothScatter_tapsbeta_oxbs_depth10.pdf")
smoothScatter(x = merged_new$oxbs_level,  y = merged_new$beta_level, xlab = "oxBS 5mC level", ylab = "TAPS-beta 5mC level")
dev.off()

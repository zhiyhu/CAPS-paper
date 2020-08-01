#!/usr/bin/env Rscript
## preprocessing CAPS data, ACE data and TABseq data
## Zhiyuan Hu
## 11 May 2020
## last modified 14 May 2020

library(data.table)
library(parallel)


## output file: ../tidy_data/5hmC_caps_mlml_ace_tab_mods.bed.gz 

## set up ----
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
} else {
  os <- "No"
}

fn_cpg  <- "$HOME/projects/taps/analysis_caps/methy_data/mm9_CpG_sites.bed.gz"
fn_caps <- "$HOME/projects/taps/analysis_caps/methy_data/caps_CpG_mods.bed.gz"
fn_ace  <- "$HOME/projects/taps/analysis_caps/methy_data/ace_CpG_mods.bed.gz"
fn_mlml <- "$HOME/projects/taps/analysis_caps/methy_data/taps_tapsbeta_mlml_CpG_mods.bed.gz"
fn_taps <- "$HOME/projects/taps/analysis_caps/methy_data/taps_CpG_mods.bed.gz"
fn_beta <- "$HOME/projects/taps/analysis_caps/methy_data/tapsbeta_CpG_mods.bed.gz"
fn_tab  <- "$HOME/projects/taps/analysis_caps/methy_data/tabseq_CpG_mods.bed.gz"
chrs <- paste("chr", c(1:19, "X","Y"), sep = "")


##################### part 1 read in files and merge -----------------

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


print("sum(is.na(cpg_all$chr)) = ")
sum(is.na(cpg_all$chr))

cpg_all <- cpg_all[!is.na(cpg_all$chr),] ## remove NA rows

fwrite(cpg_all, "../tidy_data/5hmC_caps_mlml_ace_tab_mods_binomial_tested.unfiltered.bed.gz", sep = "\t", row.names = F, col.names = T)



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



fwrite(tidy_data, "../tidy_data/5hmC_caps_mlml_ace_tab_mods.bed.gz", sep = "\t", row.names = F, col.names = T)




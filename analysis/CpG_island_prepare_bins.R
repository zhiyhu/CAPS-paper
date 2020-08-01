#!/usr/bin/env Rscript
## Analyse the modification levels and coverage at and around the CpG islands
## Zhiyuan Hu
## 20 May 2020
## last modified 1 Aug 2020

########################################################
############################  local --------------------
########################################################


## move the output files to remote projects/taps/ref/ucsc/

## input file:
## ../../data/cpgIslandExt.txt.gz

## Output files:
## cpgIslandExt_bins.bed.gz
## cpgIslandExt_left_4kb_flank_20wins.txt.gz
## cpgIslandExt_right_4kb_flank_20wins.txt.gz

### clean up the CpG island bed files -----

library(tidyr)

cgi <- fread("../../data/cpgIslandExt.txt.gz", header = F, nThread = 4,stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
colnames(cgi)[2:4] <- c("chr","start","end")
cgi <- cgi[cgi$chr %in% paste("chr", c(1:19, "X","Y"), sep = ""),]
# number of CpG islands
nrow(cgi)
# [1] 15991

cgi$length <- cgi$end - cgi$start

## cut CGIs into bins
n_bin <- 10 # number of bins
cgi$bin_width <- cgi$length/n_bin

m_bin <- matrix(0, nrow = nrow(cgi), ncol = n_bin + 1)
m_bin[,1] <- cgi$start
## make bed files for bin
for(itor in 2:(n_bin+1)){
  m_bin[,itor] <- round(m_bin[,1] + (itor - 1)*cgi$bin_width, 0)
}
m_bin <- data.frame(m_bin)
m_bin$chr <- cgi$chr
# m_bin$id <- rownames(cgi)

df1 <- gather(m_bin[,c(1:(ncol(m_bin)-2), ncol(m_bin))], bin, start, X1:X10, factor_key=TRUE)
df2 <- gather(m_bin[,c(2:(ncol(m_bin)-1), ncol(m_bin))], bin, end, X2:X11, factor_key=TRUE)
df1$end <- df2$end
df1$bin <- gsub(pattern = "X", replacement = "bin", df1$bin)
df1_bed <- data.frame(chr = df1$chr, start = df1$start, end = df1$end, name = df1$bin)
df1_bed <- df1_bed[order(df1_bed$chr, df1_bed$start),]
fwrite(df1_bed, "../../data/cpgIslandExt_bins.bed.gz", sep = "\t", col.names = F, row.names = F)

# number of bins
nrow(df1_bed)
# [1] 159910


## bin the 4kb flanks; 50 windows; 80 bp per window

## length of non-cpi regions
length <- (cgi$start[2:nrow(cgi)] - cgi$end[1:(nrow(cgi)-1)])[cgi$chr[2:nrow(cgi)] == cgi$chr[1:(nrow(cgi)-1)]]
summary(length)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 67    14184    45930   159975   127720 12190761 
sum(length<8000)
# [1] 2796

## 2796 intervals are shorter than 8kb (for 2x4kb)

# calculate the mid point of intervals (max length of the flanks)
cgi$id <- 1:nrow(cgi)
cgi$left_mid <- NA
cgi$right_mid <- NA
for(i in 2:(nrow(cgi)-1)){
  if(cgi$chr[i-1] == cgi$chr[i]) {# one the same chr with the last one
    # left middle point
    cgi$left_mid[i] <- mean(cgi$end[i-1], cgi$start[i])
  } else {
    cgi$left_mid[i] <- 0
  }
  if(cgi$chr[i] == cgi$chr[i+1]) {# one the same chr with the next one
    # right middle point
    cgi$right_mid[i] <- mean(cgi$start[i+1], cgi$end[i])
  } else {
    cgi$right_mid[i] <- cgi$end[i] + 4002
  }
}
cgi$left_mid[1] <- 0
cgi$right_mid[1] <- mean(cgi$start[2], cgi$end[1])
cgi$left_mid[nrow(cgi)] <-  mean(cgi$end[nrow(cgi)-1], cgi$start[nrow(cgi)])
cgi$right_mid[nrow(cgi)] <- cgi$end[nrow(cgi)] + 4002

## bin left flank -----
cgi$left_flank_start <- apply(cbind((cgi$start - 4002), cgi$left_mid), 1, max)

flank <-  data.frame(chr = cgi$chr,
                     start = cgi$left_flank_start,
                     end = cgi$start-1,
                     id = cgi$id)

## make bed files for bin
n_bin <- 20
length_bin <- 4000/n_bin

m_bin <- matrix(0, nrow = nrow(flank), ncol = n_bin + 1)
m_bin[,n_bin + 1] <- flank$end

for(itor in 1:(n_bin)){
  m_bin[,n_bin+1-itor] <- round(m_bin[,n_bin+1] - itor*length_bin, 0)
}
m_bin <- data.frame(m_bin)
m_bin$chr <- flank$chr
m_bin$id <- flank$id

for(i in 1:nrow(m_bin)){
  tmp <- m_bin[i, 1:21]
  m_bin[i,which(tmp < flank$start[i])] <- -5
}

df1 <- gather(m_bin[,c(1:20, 22, 23)], bin, start, X1:X20, factor_key=TRUE)
df2 <- gather(m_bin[,c(2:21, 22, 23)], bin, end, X2:X21, factor_key=TRUE)
df1$end <- df2$end
df1$bin <- gsub(pattern = "X", replacement = "bin", df1$bin)
df1_bed <- data.frame(chr = df1$chr, start = df1$start, end = df1$end, name = df1$bin, id = df1$id)
df1_bed <- df1_bed[order(df1_bed$chr, df1_bed$start),]
sum(df1_bed$start == -5)
sum(df1_bed$end == -5)

df1_bed <- df1_bed[df1_bed$start >= 0 & df1_bed$end >= 0,]

fwrite(df1_bed, "../../data/cpgIslandExt_left_4kb_flank_20wins.txt.gz", sep = "\t", col.names = F, row.names = F)

## bin right flank ----

cgi$right_flank_end <- apply(cbind((cgi$end + 4002), cgi$right_mid), 1, min)

flank <- data.frame(chr = cgi$chr,
                    start = cgi$end+1,
                    end = cgi$right_flank_end,
                    id = cgi$id)

m_bin <- matrix(0, nrow = nrow(flank), ncol = n_bin + 1)
m_bin[,1] <- flank$start+1
## make bed files for bin
for(itor in 2:ncol(m_bin)){
  m_bin[,itor] <- round(m_bin[,1] + (itor - 1)*length_bin, 0)
}
m_bin <- data.frame(m_bin)
m_bin$chr <- flank$chr
m_bin$id <- flank$id

for(i in 1:nrow(m_bin)){
  tmp <- m_bin[i, 1:21]
  m_bin[i,which(tmp > flank$end[i])] <- -5
}

df1 <- gather(m_bin[,c(1:21, 22, 23)], bin, start, X1:X20, factor_key=TRUE)
df2 <- gather(m_bin[,c(2:21, 22, 23)], bin, end, X2:X21, factor_key=TRUE)
df1$end <- df2$end
df1$bin <- gsub(pattern = "X", replacement = "bin", df1$bin)
df1_bed <- data.frame(chr = df1$chr, start = df1$start, end = df1$end, name = df1$bin, id = df1$id)
df1_bed <- df1_bed[order(df1_bed$chr, df1_bed$start),]

sum(df1_bed$start == -5)
sum(df1_bed$end == -5)

df1_bed <- df1_bed[df1_bed$start >= 0 & df1_bed$end >= 0,]

fwrite(df1_bed, "../../data/cpgIslandExt_right_4kb_flank_20wins.txt.gz", sep = "\t", col.names = F, row.names = F)
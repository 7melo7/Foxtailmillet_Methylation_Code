library(data.table)
library(data.table)
library(Biostrings)
library(pseudoRef)
args=commandArgs(T)
file1 <- args[1]
file2 <- args[2]
file3 <- args[3]
h <- read.table(file1,header=T)
snpdt <- fread(file2,header=F)
names(snpdt) <- gsub("_.*","",names(h))
fa  <- file3
pseudoRef(fa,snpdt,outdir="/public1/home/miyj/projects/sitalica/snp")



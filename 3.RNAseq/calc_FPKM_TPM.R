library(tidyverse)

args = commandArgs(T)
all_count <- read_tsv(args[1]) %>% as.data.frame()
gene_length <- all_count$Length

old_count <- read_tsv("three.counts.txt") %>% as.data.frame()
old_count <- old_count[,7:ncol(old_count)]

all_count <- select(all_count, -Chr, -Start, -End, -Strand, -Length)
rownames(all_count) <- all_count[,1]
all_count <- all_count[,-1]
all_count <- cbind(all_count, old_count)
all_count[all_count < 0] <- 0

#FPKM
length <- gene_length / 1000
fpkm <- all_count / length
fpkm <- fpkm * 1000000 / apply(all_count, 2, sum)
write_tsv(fpkm, file="all.fpkm.txt")

# TPM
tpm = all_count / length
tpm = tpm * 1000000 / apply(tpm, 2, sum)
write_tsv(tpm, file="all.tpm.txt")


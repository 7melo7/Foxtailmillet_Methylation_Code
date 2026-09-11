library(tidyverse)
library(clusterProfiler)

args <- commandArgs(T)
gene_list <- args[1]

pfam_db <- read_tsv("/public1/home/guowl/Share_GuoWL/dataset/Sitalica/func_anno/db.pfam.tsv")

set_genes <- read_tsv(gene_list, col_names = c("gene"))
set_genes <- as.vector(set_genes$gene)

longmi_pfam <- read_tsv("/public1/home/guowl/Share_GuoWL/dataset/Sitalica/func_anno/sita.pfam.tsv", col_names = c("gene", "PFAMID"))
	
pfam_term2gene <- left_join(longmi_pfam, pfam_db, by="PFAMID") %>% select("PFAMID", "gene")
pfam_term2name <- left_join(longmi_pfam, pfam_db, by="PFAMID") %>% select("PFAMID", "INFO")

names(pfam_term2gene) <- c("pfam_term", "gene")
names(pfam_term2name) <- c("pfam_term", "name")

	# pfam enrich
pfam_enrich <- enricher(gene = set_genes, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = pfam_term2gene, TERM2NAME = pfam_term2name)
pfam_enrich_result <- as_tibble(pfam_enrich@result) %>% arrange(p.adjust) %>% filter(p.adjust < 0.05, qvalue < 0.05)
write_tsv(pfam_enrich_result, file = paste0(gene_list, ".pfam.tsv"))

library(tidyverse)
library(readr)
library(clusterProfiler)

args <- commandArgs(T)
gene_list <- args[1]

# read the standard Gene Ontology and KEGG annotation database
go_db <- read_tsv("/public1/home/guowl/Share_GuoWL/dataset/standard_go_kegg/go.db.tsv")
ko_db <- read_tsv("/public1/home/guowl/Share_GuoWL/dataset/standard_go_kegg/kegg.db.tsv")

# read genes from regions
set_genes <- read_tsv(gene_list, col_names = c("gene"))
set_genes <- as.vector(set_genes$gene)
        
if(length(set_genes)==0){
    next
}
# read longmi whole genome Gene Onotology and KEGG annotation
#longmi_go <- read_tsv("~/dataset/mizi/longmi.go.tsv", col_names = c("gene", "GOID"))
longmi_go <- read_tsv("/public1/home/guowl/Share_GuoWL/dataset/Sitalica/func_anno/sita.go.interproscan.ids", col_names = c("gene", "GOID"))
#longmi_go <- read_tsv("~/dataset/mizi/longmi.go.tsv", col_names = c("gene", "GOID"))
longmi_ko <- read_tsv("/public1/home/guowl/Share_GuoWL/dataset/Sitalica/func_anno/sita.ko.hmm.ids", col_names = c("gene", "KO"))

# prepare enricher data
go_term2gene <- left_join(longmi_go, go_db, by="GOID") %>% select("GOID", "gene")
go_term2name <- left_join(longmi_go, go_db, by="GOID") %>% select("GOID", "TERM")

names(go_term2gene) <- c("go_term", "gene")
names(go_term2name) <- c("go_term", "name")

kegg_term2gene <- left_join(longmi_ko, ko_db, by="KO") %>% select("KO", "gene")
kegg_term2name <- left_join(longmi_ko, ko_db, by="KO") %>% select("KO", "INFO")

names(kegg_term2gene) <- c("ko_term", "gene")
names(kegg_term2name) <- c("ko_term", "name")

# go enrich
go_enrich <- enricher(gene = set_genes, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     TERM2GENE = go_term2gene, TERM2NAME = go_term2name)
if(is.null(go_enrich)){
    next
}
go_enrich_result <- as_tibble(go_enrich@result) %>% arrange(p.adjust) %>% filter(p.adjust < 0.05, qvalue < 0.05)

write_tsv(go_enrich_result, file = paste0(gene_list, ".go.enrich.tsv"))

#pdf(file=("go.pdf"), width=15, height=15)
#emapplot(go_enrich, showCategory=50)
#dev.off()
# go DAG plot
#pdf(file=paste0("enrichment/plot/", j, ".win_", i, ".go.bp.tree.pdf"), width=10, height=15)
#plotGOgraph(go_enrich)
#dev.off()

#pdf(file=paste0("enrichment/plot/", j, ".win_", i, ".go.mf.tree.pdf"), width=10, height=15)
#plotGOgraph(go_enrich.MF)
#dev.off()

#pdf(file=paste0("enrichment/plot/", j, ".win_", i, ".go.cc.tree.pdf"), width=10, height=15)
#plotGOgraph(go_enrich.CC)
#dev.off()
#barplot(c3_go_enrich, showCategory = 40)

# kegg enrich
ko_enrich <- enricher(gene = set_genes, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name)
if (is.null(ko_enrich)){
    next
}
#pdf(file=("ko.pdf"), width=15, height=15)
#emapplot(ko_enrich, showCategory=50)
#dev.off()
ko_enrich_result <- as_tibble(ko_enrich@result) %>% arrange(p.adjust) %>% filter(p.adjust < 0.05, qvalue < 0.05)

write_tsv(ko_enrich_result, file = paste0(gene_list, ".ko.enrich.tsv"))
#barplot(all_ko_enrich, showCategory = 40)






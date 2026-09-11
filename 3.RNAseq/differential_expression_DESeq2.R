library("DESeq2")
library("tidyverse")
library("ggplot2")

args <- commandArgs(T)

all_count <- read.table(args[1], header = T, row.names = 1)
all_count <- as.matrix(all_count)
all_count <- all_count[rowSums(all_count) > ncol(all_count),]


col_group <- tibble(samples = colnames(all_count))
col_group$pops <- sapply(col_group$samples, function(x){
    if(str_starts(x, 'Q')){
        return("wild")
    }else{
        return("cultivar")
    }
})
col_group$pops <- as.factor(col_group$pops)

dds <- DESeqDataSetFromMatrix(countData = all_count,
                              colData = col_group,
                              design = ~ pops)
# PCA plot

vsd <- vst(dds, blind = FALSE)
vsd_pca <- plotPCA(vsd, intgroup = c("pops"), returnData = T)

library("ggpubr")
library("ggthemes")

p <- ggscatter(vsd_pca, x="PC1", y="PC2", shape="pops", color="pops",
               palette = c("#548235", "#4472C4"), size=2, ellipse.border.remove=FALSE, ellipse=TRUE,
               repel=TRUE, main="PCA plot"+theme_base())

pdf("rnaseq_pca.pdf", width = 5, height = 5)
p
dev.off()


dds <- DESeq(dds)
resultsNames(dds)
dds_res <- results(dds, name="pops_wild_vs_cultivar")
dds_res <- as.tibble(dds_res)
dds_res$gene <- row.names(all_count)
dds_res <- select(dds_res, gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
diff_dds_res <- filter(dds_res, padj < 0.05, abs(log2FoldChange) > 1)

write_tsv(diff_dds_res, "pop_diffed.gene_info.tsv")
write_tsv(dds_res, "pop_diff.gene_info.tsv")

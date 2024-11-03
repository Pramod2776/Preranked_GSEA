setwd("~/Path/to/Workdir/")
library(tidyverse)
res = read_tsv("./deseq2.results.tsv")

###  Add annotation to deseq2 results
geneanno <- read.delim("./geneInfo.tab", header = F,
                       stringsAsFactors = F, skip =1)
geneanno <- geneanno[,1:3]
colnames(geneanno) <- c("gene_id", "gene_name", "gene_type")

# make gene names unique
genenames <- geneanno$gene_name
genenames[which(duplicated(genenames))]
which(duplicated(genenames))
genenames.fixed <- genenames
genenames.fixed[which(duplicated(genenames.fixed))] <- paste0(genenames.fixed[which(duplicated(genenames.fixed))],"_",
                                                              geneanno$gene_id[which(duplicated(genenames.fixed))])
genenames.fixed[which(duplicated(genenames))]
geneanno$gene_name_unique <- genenames.fixed

res_anno = res %>%
  dplyr::left_join(geneanno, by = "gene_id") %>%
  dplyr::select("gene_id", "gene_type", "gene_name", "gene_name_unique",
                "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj") %>%
  as.data.frame()
res_anno[res_anno ==""]<-NA

### make ranked .rnk file for preranked gsea
## In the following code we generate a ranked data frame with the gene
## names as the first column. The second column is the negative log-p
## value multiplied by the sign of the fold change.
rnkdf <- tibble(gene = res_anno$gene_name_unique,
                rnk = -log(res_anno$pvalue) * sign(res_anno$log2FoldChange)) %>%
  arrange(desc(rnk)) %>% drop_na()

## Write out the table without any additional information
write.table(rnkdf, file = "deseq2_res_for_gsea_outlier.rnk",append = FALSE, col.names = FALSE, row.names = FALSE,
quote = FALSE, sep = "\t")

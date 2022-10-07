
# read data
setwd('/data/OriginalData/TCGA-LIHC')
tcga.lihc.exp <- read.csv(file='LIHC_rnaseqv2_RSEM_genes.txt', sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE)

tcga.lihc.count <- tcga.lihc.exp[-1, seq(1, ncol(tcga.lihc.exp), 3)]
tcga.lihc.tpm <- tcga.lihc.exp[-1, seq(2, ncol(tcga.lihc.exp), 3)]

colnames(tcga.lihc.count) <- substr(gsub('\\.', '-', colnames(tcga.lihc.count)), 1, 15)
colnames(tcga.lihc.tpm) <- substr(gsub('\\.', '-', colnames(tcga.lihc.tpm)), 1, 15)


# Raw count and scaled estimate (TPM) (PMID:34430923)
# We extracted the gene expression data from "illuminahiseq_rnaseqv2-RSEM_genes" files.
# From these data, we used "raw_count" values as counts and 
# we calculated transcripts per million (TPM) from "scaled_estimate" values multiplied by 1,000,000. 


library(dplyr)
tmp <- tcga.lihc.count %>% mutate(across(where(is.character), as.numeric))
rownames(tmp) <- rownames(tcga.lihc.count)
tcga.lihc.count <- tmp
tcga.lihc.count <- round(tcga.lihc.count)

tmp <- tcga.lihc.tpm %>% mutate(across(where(is.character), as.numeric))
rownames(tmp) <- rownames(tcga.lihc.tpm)
tcga.lihc.tpm <- tmp
tcga.lihc.tpm <- tcga.lihc.tpm * 10^6


save(tcga.lihc.tpm, tcga.lihc.count, file='/data/tcga_lihc_rnaseq.RData')


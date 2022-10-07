
setwd('/data/OriginalData/GSE144269')
cli.data <- read.csv(file='patient_sample_metadata.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
cli.data <- subset(cli.data, !is.na(RNASeq_T))


exp.count <- read.csv(file='GSE144269_RSEM_GeneCounts.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)

exp.count$entrez_id <- do.call(rbind, strsplit(exp.count$entrez_id, split='\\|'))[, 1]
exp.count <- tibble::column_to_rownames(exp.count, var = "entrez_id")

exp.count <- exp.count[rowSums(exp.count >= 1) >= 28, ] 


# sample information
sams.info <- data.frame(labels=c(rep('Tumor', nrow(cli.data)), rep('Normal', nrow(cli.data))))
rownames(sams.info) <- c(cli.data$RNASeq_T, cli.data$RNASeq_NT)

exp.count <- exp.count[, rownames(sams.info)]

# Differential expression analysis
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData = sams.info, design = ~ labels)

vsd <- vst(dds, blind = FALSE)
GSE144269.vst.exp <- assay(vsd)

dds <- DESeq(dds)
GSE144269.diff.res <- results(dds)

save(cli.data, sams.info, GSE144269.vst.exp, GSE144269.diff.res, file='/data/GSE144269_diff_exp.RData')


# candi.prog.gene <- c("ENSG00000116128", "ENSG00000129993", "ENSG00000147883", "ENSG00000102967",
# "ENSG00000187486", "ENSG00000104738", "ENSG00000175806", "ENSG00000198056",
# "ENSG00000125249", "ENSG00000197747", "ENSG00000137642", "ENSG00000234616",
# "ENSG00000164741", "ENSG00000160767", "ENSG00000184792", "ENSG00000138744",
# "ENSG00000161800", "ENSG00000146918", "ENSG00000101695", "ENSG00000150403",
# "ENSG00000137269", "ENSG00000025039", "ENSG00000225697", "ENSG00000160050",
# "ENSG00000185803", "ENSG00000120262", "ENSG00000162174", "ENSG00000157600",
# "ENSG00000181751", "ENSG00000158865", "ENSG00000149639", "ENSG00000139914",
# "ENSG00000165633", "ENSG00000188368")

# > as.data.frame(GSE144269.diff.res[pmatch(candi.prog.gene, rownames(GSE144269.diff.res)),])$padj
 # [1] 5.104160e-25 7.019193e-07 2.730263e-36 3.652990e-20 5.225280e-03 4.684027e-34 8.284740e-27 1.659189e-21 5.310626e-22
# [10] 1.190514e-22 1.066757e-23 4.126338e-08 7.581788e-06 3.988829e-31 3.836873e-13 6.710544e-17 1.841339e-50 2.164758e-36
# [19] 1.482477e-12 1.245613e-20 1.124680e-13 3.418282e-27 1.662900e-26 3.051380e-10 1.066105e-13 1.413073e-18 1.624469e-17
# [28] 7.971496e-12 5.712540e-12 1.900310e-07 8.785962e-20 8.071464e-16 1.636963e-02 2.709685e-08



# setwd('/result/Section2/')
# colnames(sams.info) <- 'Type'
# GSE144269.vst.exp <- GSE144269.vst.exp[!duplicated(substr(rownames(GSE144269.vst.exp), 1, 15)), ]
# rownames(GSE144269.vst.exp) <- substr(rownames(GSE144269.vst.exp), 1, 15)
# ExpClusteringHeatmap(GSE144269.vst.exp, candi.prog.gene, sams.info, 'GSE144269_candi_gene_exp_heatmap')




# load data
setwd('/data/OriginalData/GSE144269')
cli.data <- read.csv(file='patient_sample_metadata.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
cli.data <- subset(cli.data, !is.na(RNASeq_T))


exp.count <- read.csv(file='GSE144269_RSEM_GeneCounts.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)

exp.count$entrez_id <- do.call(rbind, strsplit(exp.count$entrez_id, split='\\|'))[, 1]
exp.count <- tibble::column_to_rownames(exp.count, var = "entrez_id")

exp.count <- exp.count[, cli.data$RNASeq_T]
exp.count <- exp.count[rowSums(exp.count >= 1) >= 14, ] 


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=cli.data, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
norm.exp <- assay(vsd)

# save
saveRDS(cli.data, file='/data/GSE144269_expr_data.rds')
saveRDS(norm.exp, file='/data/GSE144269_cli_data.rds')

# load data
setwd('/data/OriginalData/Gao_Cell_2019')
identifier <- readxl::read_xlsx(path='About the RNA and protein Identifier match.xlsx', sheet = 1, col_names = TRUE, skip=1)
colnames(identifier) <- c('protein_T', 'protein_N', 'rna_T', 'rna_N')


# expression data
exp.data <- read.csv(file='HCC_UQ_RSEM.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
exp.data <- tibble::column_to_rownames(exp.data, var = "protein")
exp.data <- dplyr::mutate(exp.data, X=NULL)
colnames(exp.data) <- gsub('X', '', colnames(exp.data))

exp.data <- exp.data[, as.character(c(identifier$rna_T, identifier$rna_N))]
colnames(exp.data) <- c(paste0('T', identifier$rna_T), paste0('N', identifier$rna_N))


exp.data <- round(exp.data)
exp.data <- exp.data[rowSums(exp.data >= 1) >= 64, ]



# sample information
sams.info <- data.frame(labels=c(rep('Tumor', nrow(identifier)), rep('Normal', nrow(identifier))))
rownames(sams.info) <- c(paste0('T', identifier$rna_T), paste0('N', identifier$rna_N))


# Differential expression analysis
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exp.data, colData = sams.info, design = ~ labels)

vsd <- vst(dds, blind = FALSE)
gao.etal.vst.exp <- assay(vsd)

dds <- DESeq(dds)
gao.etal.diff.res <- results(dds)

save(sams.info, gao.etal.vst.exp, gao.etal.diff.res, file='/data/Gao_Cell_2019_diff_exp.RData')


# candi.prog.gene <- c("ENSG00000116128", "ENSG00000129993", "ENSG00000147883", "ENSG00000102967",
# "ENSG00000187486", "ENSG00000104738", "ENSG00000175806", "ENSG00000198056",
# "ENSG00000125249", "ENSG00000197747", "ENSG00000137642", "ENSG00000234616",
# "ENSG00000164741", "ENSG00000160767", "ENSG00000184792", "ENSG00000138744",
# "ENSG00000161800", "ENSG00000146918", "ENSG00000101695", "ENSG00000150403",
# "ENSG00000137269", "ENSG00000025039", "ENSG00000225697", "ENSG00000160050",
# "ENSG00000185803", "ENSG00000120262", "ENSG00000162174", "ENSG00000157600",
# "ENSG00000181751", "ENSG00000158865", "ENSG00000149639", "ENSG00000139914",
# "ENSG00000165633", "ENSG00000188368")

# as.data.frame(gao.etal.diff.res[pmatch(candi.prog.gene, rownames(gao.etal.diff.res)),])$padj
# 1.278957e-69 1.295259e-104  1.221313e-38  4.397604e-29  1.822098e-17 1.156417e-105  6.090486e-47  5.352033e-81  1.816262e-68
# 1.329217e-48  3.778584e-46  2.071378e-28  4.146505e-24  7.171724e-96  3.304486e-44  2.525463e-92 8.435374e-154 5.978125e-104
# 7.066051e-38  6.788744e-59  7.021549e-60  5.827771e-72  2.632762e-42  2.618960e-35  1.544743e-36  1.522603e-44  3.572155e-70
# 7.652955e-46  1.399946e-32  3.874562e-21  4.750621e-58  1.416953e-58  7.738024e-01  8.230638e-17


# setwd('/result/Section2/')
# colnames(sams.info) <- 'Type'
# gao.etal.vst.exp <- gao.etal.vst.exp[!duplicated(substr(rownames(gao.etal.vst.exp), 1, 15)), ]
# rownames(gao.etal.vst.exp) <- substr(rownames(gao.etal.vst.exp), 1, 15)
# ExpClusteringHeatmap(gao.etal.vst.exp, candi.prog.gene, sams.info, 'gao.etal_candi_gene_exp_heatmap')


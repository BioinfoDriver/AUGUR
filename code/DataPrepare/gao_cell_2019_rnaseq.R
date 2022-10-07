

# load data
setwd('/data/OriginalData/Gao_Cell_2019')
cli.data <- readxl::read_xlsx(path='mmc1.xlsx', sheet = 2, col_names = TRUE)
cli.data <- as.data.frame(cli.data)
rownames(cli.data) <- cli.data$`Tumor (T) sample ID`

cli.data$`Overall survial (month)` <- cli.data$`Overall survial (month)`*30.4375
cli.data$`Recurrence-free survival (month)` <- cli.data$`Recurrence-free survival (month)`*30.4375



identifier <- readxl::read_xlsx(path='About the RNA and protein Identifier match.xlsx', sheet = 1, col_names = TRUE, skip=1)
colnames(identifier) <- c('protein_T', 'protein_N', 'rna_T', 'rna_N')

# expression data
exp.data <- read.csv(file='HCC_UQ_RSEM.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
exp.data <- tibble::column_to_rownames(exp.data, var = "protein")
exp.data <- dplyr::mutate(exp.data, X=NULL)
colnames(exp.data) <- gsub('X', '', colnames(exp.data))


exp.t.data <- exp.data[, as.character(identifier$rna_T)]
colnames(exp.t.data) <- paste0('T', colnames(exp.t.data))

exp.data <- round(exp.t.data)
exp.data <- exp.data[rowSums(exp.data >= 1) >= 32, ]



library(DESeq2)
col.data <- data.frame(label=c(rep('Tumor', nrow(cli.data))))
rownames(col.data) <- c(cli.data$`Tumor (T) sample ID`)

dds <- DESeqDataSetFromMatrix(countData = exp.data, colData=col.data, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
norm.exp <- assay(vsd)

# save
saveRDS(cli.data, file='/data/Gao_Cell_2019_expr_data.rds')
saveRDS(norm.exp, file='/data/Gao_Cell_2019_cli_data.rds')
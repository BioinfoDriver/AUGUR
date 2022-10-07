
# load data
tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')
tcga.lihc.cli.data <- readRDS(file='/data/tcga.lihc.cli.data.rds')
load('/data/exp_gene_anno.RData')

tcga.lihc.vst.exp <- tcga.lihc.vst.exp[as.character(na.omit(exp.gene.anno$TCGAExpGeneFilter)), ]

# data preparation
rownames(tcga.lihc.vst.exp) <- paste0('ID_', rownames(tcga.lihc.vst.exp))
colnames(tcga.lihc.vst.exp) <- substr(colnames(tcga.lihc.vst.exp), 1, 12)
tcga.lihc.vst.exp <- as.data.frame(t(tcga.lihc.vst.exp))


tcga.lihc.cli.data <- tcga.lihc.cli.data[, c('os_time', 'os')]
cli.exp.data <- merge(tcga.lihc.cli.data, tcga.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})

names(univ.cox.pvalue) <- gsub('ID_', '', names(univ.cox.pvalue))
saveRDS(univ.cox.pvalue, file='/data/tcga.lihc.univ.cox.pvalue.rds')

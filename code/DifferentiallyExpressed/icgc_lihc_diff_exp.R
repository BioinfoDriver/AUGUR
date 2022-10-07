
# load data
icgc.rnaseq.sample <- readRDS('/data/icgc_rnaseq_filter_sample.rds')
load(file='/data/icgc_lihc_filter_rnaseq.RData')
load(file='/data/exp_gene_anno.RData')

icgc.lihc.filter.count <- icgc.lihc.filter.count[na.omit(exp.gene.anno$ICGCExpGeneFilter), ]
rownames(icgc.lihc.filter.count) <- exp.gene.anno$GeneID[match(rownames(icgc.lihc.filter.count), exp.gene.anno$ICGCExpGeneFilter)]

#replace all NA values with zero
library(dplyr)
icgc.lihc.filter.count <- icgc.lihc.filter.count %>% replace(is.na(.), 0)


# sample information
icgc.rnaseq.sample$submitted_sample_id <- do.call(rbind, strsplit(icgc.rnaseq.sample$submitted_sample_id, split='_'))[, 1]
tumor.sams <- subset(icgc.rnaseq.sample, specimen_type == 'Tumor')$submitted_sample_id
normal.sams <- subset(icgc.rnaseq.sample, specimen_type == 'Normal')$submitted_sample_id
paired.sams <- setdiff(intersect(tumor.sams, normal.sams), c('RK062', 'RK080')) # 173


# primary and normal samples
tumor.sams <- subset(icgc.rnaseq.sample, submitted_sample_id %in% paired.sams & specimen_type == 'Tumor')$icgc_sample_id
normal.sams <- subset(icgc.rnaseq.sample, submitted_sample_id %in% paired.sams & specimen_type == 'Normal')$icgc_sample_id
paired.sams.count <- icgc.lihc.filter.count[, c(tumor.sams, normal.sams)]


# sample information
paired.sams.info <- data.frame(labels=c(rep('Tumor', length(paired.sams)), rep('Normal', length(paired.sams))))
rownames(paired.sams.info) <- c(tumor.sams, normal.sams)


# Differential expression analysis
library(DESeq2)
paired.dds <- DESeqDataSetFromMatrix(countData = paired.sams.count, colData = paired.sams.info, design = ~ labels)
paired.dds <- DESeq(paired.dds)
icgc.paired.diff.res <- results(paired.dds)


# save
save(icgc.paired.diff.res, file='/data/icgc_lihc_diff_exp.RData')


# load data
load(file='/data/tcga_lihc_rnaseq.RData')
load(file='/data/exp_gene_anno.RData')
load(file='/data/tcga_lihc_tv_samples.RData')


rownames(tcga.lihc.count) <- do.call(rbind, strsplit(rownames(tcga.lihc.count), split='\\|'))[, 2]
tcga.lihc.count <- tcga.lihc.count[as.character(na.omit(exp.gene.anno$TCGAExpGeneFilter)), ]


# primary and normal samples 
normal.sams <- colnames(tcga.lihc.count)[substr(colnames(tcga.lihc.count), 14, 15) == '11'] # 50
normal.sams <- substr(normal.sams, 1, 12)
primary.sams <- colnames(tcga.lihc.count)[substr(colnames(tcga.lihc.count), 14, 15) == '01'] # 371
primary.sams <- substr(primary.sams, 1, 12)
paired.sams <- intersect(normal.sams, primary.sams) # 50

paired.sams <- intersect(paired.sams, c(tcga.train.sam.set, tcga.test.sam.set)) # 42
primary.sams <- intersect(primary.sams, c(tcga.train.sam.set, tcga.test.sam.set)) # 323

paired.sams.count <- tcga.lihc.count[, c(paste0(paired.sams, '-01'), paste0(paired.sams, '-11'))]
unpaired.sams.count <- tcga.lihc.count[, c(paste0(primary.sams, '-01'), paste0(paired.sams, '-11'))] 

# sample information
paired.sams.info <- data.frame(labels=c(rep('Tumor', length(paired.sams)), rep('Normal', length(paired.sams))))
rownames(paired.sams.info) <- c(paste0(paired.sams, '-01'), paste0(paired.sams, '-11'))

unpaired.sams.info <- data.frame(labels=c(rep('Tumor', length(primary.sams)), rep('Normal', length(paired.sams))))
rownames(unpaired.sams.info) <- c(primary.sams, paired.sams)


# Differential expression analysis
library(DESeq2)
paired.dds <- DESeqDataSetFromMatrix(countData = paired.sams.count, colData = paired.sams.info, design = ~ labels)
unpaired.dds <- DESeqDataSetFromMatrix(countData = unpaired.sams.count, colData = unpaired.sams.info, design = ~ labels)
# save(paired.dds, unpaired.dds, file='/data/tcga_lihc_exp_dds.RData')


paired.dds <- DESeq(paired.dds)
unpaired.dds <- DESeq(unpaired.dds)

tcga.paired.diff.res <- results(paired.dds)
tcga.unpaired.diff.res <- results(unpaired.dds)


# save data
save(tcga.paired.diff.res, tcga.unpaired.diff.res, file='/data/tcga_lihc_diff_exp.RData')


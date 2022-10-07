
# load data
# Gene expression from Losic et al.
load(file='/data/losic.refseq.gene.exp.RData')
load(file='/data/exp_gene_anno.RData')

# sample information
setwd('/data/OriginalData/Losic_Nat_Commun_2020_LiverMultiregionExp')
sam.info <- read.csv(file='E-MTAB-5905.sdrf.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
rownames(sam.info) <- sam.info$Source.Name
sam.info <- subset(sam.info, Characteristics.sampling.site. == 'neoplasm')

# keeping genes with an expression value of at least 1 TPM in at least 20% (9/43) of tumor samples
losic.exp.tpm <- tran.exp$abundance[, rownames(sam.info)]
rownames(losic.exp.tpm) <- do.call(rbind, strsplit(rownames(losic.exp.tpm), split='_'))[, 2]
losic.exp.genes <- rownames(losic.exp.tpm)[rowSums(losic.exp.tpm >= 1) >= 9]

# length(losic.exp.genes)/nrow(losic.exp.tpm)——14401/19317=0.746
# length(intersect(losic.exp.genes, na.omit(exp.gene.anno$LosicExpGene)))/length(na.omit(exp.gene.anno$LosicExpGene))——14321/19170=0.747

losic.exp.count <- tran.exp$counts[, rownames(sam.info)]
rownames(losic.exp.count) <- do.call(rbind, strsplit(rownames(losic.exp.count), split='_'))[, 2]

# length(rownames(losic.exp.count)[rowSums(losic.exp.count >= 1) >= 9]), 17717
# length(rownames(losic.exp.count)[rowSums(losic.exp.count >= 10) >= 9]), 15663
# length(rownames(losic.exp.count)[rowSums(losic.exp.count >= 100) >= 9]), 12721



# TCGA LIHC's gene expression
load(file='/data/tcga_lihc_rnaseq.RData')
load(file='/data/tcga_lihc_tv_samples.RData')

# keeping genes with an expression value of at least 1 TPM in at least 20% (65/323) of tumor samples
rownames(tcga.lihc.tpm) <- sapply(strsplit(rownames(tcga.lihc.tpm), split='\\|'), function(x) x[2])
tcga.lihc.tpm <- tcga.lihc.tpm[, paste0(c(tcga.train.sam.set, tcga.test.sam.set), '-01')]
tcga.exp.genes <- rownames(tcga.lihc.tpm)[rowSums(tcga.lihc.tpm >= 1) >= 65] 

# length(tcga.exp.genes)/nrow(tcga.lihc.tpm)——13432/20531=0.654
# length(intersect(tcga.exp.genes, na.omit(exp.gene.anno$TCGAExpGene)))/length(na.omit(exp.gene.anno$TCGAExpGene))——12895/18534=0.696

# rownames(tcga.lihc.count) <- sapply(strsplit(rownames(tcga.lihc.count), split='\\|'), function(x) x[2])
# tcga.lihc.count <- tcga.lihc.count[, paste0(c(tcga.train.sam.set, tcga.test.sam.set), '-01')]

# length(rownames(tcga.lihc.count)[rowSums(tcga.lihc.count >= 1) >= 65]), 18204
# length(rownames(tcga.lihc.count)[rowSums(tcga.lihc.count >= 10) >= 65]), 15875
# length(rownames(tcga.lihc.count)[rowSums(tcga.lihc.count >= 100) >= 65]), 13139



# ICGC LIHC's gene expression
load(file='/data/icgc_lihc_filter_rnaseq.RData')
sam.info <- readRDS('/data/icgc_rnaseq_filter_sample.rds')
sam.info <- subset(sam.info, specimen_type == 'Tumor')

library(dplyr)
icgc.lihc.filter.fpkm <- icgc.lihc.filter.fpkm %>% replace(is.na(.), 0)
icgc.lihc.filter.count <- icgc.lihc.filter.count %>% replace(is.na(.), 0)
icgc.lihc.filter.fpkm <- icgc.lihc.filter.fpkm[, sam.info$icgc_sample_id]
icgc.lihc.filter.count <- icgc.lihc.filter.count[, sam.info$icgc_sample_id]

# keeping genes with an expression value of at least 1 read in at least 20% (41/203) of tumor samples
# icgc.exp.genes <- rownames(icgc.lihc.filter.fpkm)[rowSums(icgc.lihc.filter.fpkm >= 1) >= 41] 

# length(icgc.exp.genes)/nrow(icgc.lihc.filter.fpkm)——14111/22913=0.616
# length(intersect(icgc.exp.genes, na.omit(exp.gene.anno$ICGCExpGene)))/length(na.omit(exp.gene.anno$ICGCExpGene))——12958/18634=0.695

icgc.exp.genes <- rownames(icgc.lihc.filter.count)[rowSums(icgc.lihc.filter.count >= 1) >= 41] 

# length(icgc.exp.genes)/nrow(icgc.lihc.filter.count)——20276/22913=0.0.885
# length(intersect(icgc.exp.genes, na.omit(exp.gene.anno$ICGCExpGene)))/length(na.omit(exp.gene.anno$ICGCExpGene))——17289/18634=0.928


# length(rownames(icgc.lihc.filter.count)[rowSums(icgc.lihc.filter.count >= 1) >= 41]), 20276
# length(rownames(icgc.lihc.filter.count)[rowSums(icgc.lihc.filter.count >= 10) >= 41]), 20237
# length(rownames(icgc.lihc.filter.count)[rowSums(icgc.lihc.filter.count >= 100) >= 41]), 19987
# length(rownames(icgc.lihc.filter.count)[rowSums(icgc.lihc.filter.count >= 1000) >= 41]), 17730
# length(rownames(icgc.lihc.filter.count)[rowSums(icgc.lihc.filter.count >= 10000) >= 41]), 14724

# save
exp.gene.anno$LosicExpGeneFilter <- exp.gene.anno$LosicExpGene
exp.gene.anno$TCGAExpGeneFilter <- exp.gene.anno$TCGAExpGene
exp.gene.anno$ICGCExpGeneFilter <- exp.gene.anno$ICGCExpGene

exp.gene.anno$LosicExpGeneFilter[!(exp.gene.anno$LosicExpGeneFilter %in% losic.exp.genes)] <- NA #14321
exp.gene.anno$TCGAExpGeneFilter[!(exp.gene.anno$TCGAExpGeneFilter %in% tcga.exp.genes)] <- NA # 12895
exp.gene.anno$ICGCExpGeneFilter[!(exp.gene.anno$ICGCExpGeneFilter %in% icgc.exp.genes)] <- NA # 17289

save(exp.gene.anno, file='/data/exp_gene_anno.RData') # 12702





# load data
load(file='/data/tcga_lihc_rnaseq.RData')
load(file='/data/icgc_lihc_rnaseq.RData')
load(file='/data/losic.refseq.gene.exp.RData')

load(file='/data/exp_gene_anno.RData')

library(DESeq2)
# The variance stabilizing transformation——TCGA
load(file='/data/tcga_lihc_tv_samples.RData')
rownames(tcga.lihc.count) <- do.call(rbind, strsplit(rownames(tcga.lihc.count), split='\\|'))[, 2]
tcga.lihc.count <- tcga.lihc.count[, paste0(c(tcga.train.sam.set, tcga.test.sam.set), '-01')]
tcga.lihc.count <- tcga.lihc.count[as.character(na.omit(exp.gene.anno$TCGAExpGene)), ]

# sample information
sam.info <- data.frame(sample.type=substr(colnames(tcga.lihc.count), 14, 15))
rownames(sam.info) <- colnames(tcga.lihc.count)

dds <- DESeqDataSetFromMatrix(countData = tcga.lihc.count, colData=sam.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
tcga.lihc.vst.exp <- assay(vsd)



# The variance stabilizing transformation——ICGC
icgc.lihc.count <- icgc.lihc.count[na.omit(exp.gene.anno$ICGCExpGene), ]
rownames(icgc.lihc.count) <- exp.gene.anno$GeneID[match(rownames(icgc.lihc.count), exp.gene.anno$ICGCExpGene)]

icgc.rnaseq.sample <- readRDS(file='/data/icgc_rnaseq_filter_sample.rds')
icgc.rnaseq.sample <- subset(icgc.rnaseq.sample, specimen_type == 'Tumor')
icgc.lihc.count <- icgc.lihc.count[, icgc.rnaseq.sample$icgc_sample_id]

sam.info <- data.frame(sample.type=icgc.rnaseq.sample$specimen_type)
rownames(sam.info) <- icgc.rnaseq.sample$icgc_sample_id

icgc.lihc.count <- icgc.lihc.count %>% replace(is.na(.), 0)

dds <- DESeqDataSetFromMatrix(countData = icgc.lihc.count, colData=sam.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
icgc.lihc.vst.exp <- assay(vsd)



# The variance stabilizing transformation——Losic
# sample information
setwd('/data/OriginalData/Losic_Nat_Commun_2020_LiverMultiregionExp')
sam.info <- read.csv(file='E-MTAB-5905.sdrf.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

rownames(sam.info) <- sam.info$Source.Name
sam.info <- subset(sam.info, Characteristics.sampling.site. == 'neoplasm')

losic.exp.count <- tran.exp$counts[, rownames(sam.info)]

# duplicated
rownames(losic.exp.count) <- do.call(rbind, strsplit(rownames(losic.exp.count), split='_'))[, 2]
losic.exp.count <- losic.exp.count[-(which(duplicated(rownames(losic.exp.count)))), ]
losic.exp.count <- losic.exp.count[na.omit(exp.gene.anno$LosicExpGene), ]
rownames(losic.exp.count) <- exp.gene.anno$GeneID[match(rownames(losic.exp.count), exp.gene.anno$LosicExpGene)]


dds <- DESeqDataSetFromMatrix(countData = round(losic.exp.count), colData=sam.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
losic.lihc.vst.exp <- assay(vsd)


saveRDS(tcga.lihc.vst.exp, file='/data/tcga_vst_norm_geneExp.rds')
saveRDS(icgc.lihc.vst.exp, file='/data/icgc_vst_norm_geneExp.rds')
saveRDS(losic.lihc.vst.exp, file='/data/losic_vst_norm_geneExp.rds')

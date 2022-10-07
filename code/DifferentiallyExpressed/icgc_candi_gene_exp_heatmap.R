
load(file='/data/exp_gene_anno.RData') 
load(file='/data/icgc_lihc_diff_exp.RData')
candi.prog.gene <- readRDS(file='/data/candi_prog_sig_gene.rds')
icgc.diff.gene <- as.data.frame(icgc.paired.diff.res[candi.prog.gene, ])

# diff gene
signature <- exp.gene.anno$Symbol[match(candi.prog.gene, exp.gene.anno$GeneID)]
rownames(icgc.diff.gene) <- signature



# heatmap
# load data
icgc.rnaseq.sample <- readRDS('/data/icgc_rnaseq_filter_sample.rds')
load(file='/data/icgc_lihc_filter_rnaseq.RData')

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


library(DESeq2)
paired.dds <- DESeqDataSetFromMatrix(countData = paired.sams.count, colData = paired.sams.info, design = ~ labels)
vsd <- vst(paired.dds, blind = FALSE)
vst.exp <- assay(vsd)
rownames(vst.exp) <- exp.gene.anno$Symbol[match(rownames(vst.exp), exp.gene.anno$GeneID)]



ExpClusteringHeatmap <- function(exp.data, gene.sig, sample.anno, file.name){
  library(ComplexHeatmap)
  library(dplyr)
  library(tibble)
  library(circlize)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)
  

  # mat1
  mat1 <- exp.data[which(rownames(exp.data) %in% gene.sig), ] %>% t() %>% as.data.frame() 

  # Z-score
  tmp <- sapply(X=mat1, FUN=scale)
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- rownames(mat1)
  mat1 <- tmp
  
  
  # Annotation
  sample.anno <- sample.anno[rownames(mat1), ]
  row.col = structure(names = c('Tumor', 'Normal'), c("black", "blue"))
  row.ha <- HeatmapAnnotation(df=sample.anno, name = 'Type', which='row', col=list(Type=row.col))
  

  # ht1: expression heatmap
  # make heatmap
  ht1 = Heatmap(
    matrix=as.matrix(mat1), 
    name = "ht1", 
    col = colorRamp2(seq(-3,3), viridis(7, option = "magma")),
    column_dend_height = unit(15, "mm"),
    row_dend_width = unit(15, "mm"),
    show_row_names = FALSE,
    width = 1,
    heatmap_legend_param = list(title = NULL, color_bar = "continuous"),
    column_names_gp = gpar(fontsize = 6), 
    cluster_columns = TRUE,
	left_annotation = row.ha
  )
  
  
  # plot ComplexHeatmaps 
  pdf(paste0(file.name,'.pdf'))
  print(ht1)
  
  # ht1: add border
  decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})

  dev.off()
}


setwd('/result/Section2')

sam.anno <- paired.sams.info
colnames(sam.anno) <- 'Type'
ExpClusteringHeatmap(vst.exp, signature, sam.anno, 'icgc_candi_gene_exp_heatmap')


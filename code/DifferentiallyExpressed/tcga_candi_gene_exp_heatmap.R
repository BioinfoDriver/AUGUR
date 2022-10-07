

load(file='/data/exp_gene_anno.RData') 

candi.prog.gene <- readRDS(file='/data/candi_prog_sig_gene.rds')
signature <- exp.gene.anno$Symbol[match(candi.prog.gene, exp.gene.anno$TCGAExpGene)]


load(file='/data/tcga_lihc_exp_dds.RData')
vsd <- vst(paired.dds, blind = FALSE)
vst.exp <- assay(vsd)
rownames(vst.exp) <- exp.gene.anno$Symbol[match(rownames(vst.exp), exp.gene.anno$TCGAExpGene)]

sam.anno <- data.frame(Type=ifelse(substr(colnames(vst.exp), 14, 15) == '01',  'Tumor', 'Normal'))
rownames(sam.anno) <- colnames(vst.exp)


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
ExpClusteringHeatmap(vst.exp, signature, sam.anno, 'tcga_candi_gene_exp_heatmap')


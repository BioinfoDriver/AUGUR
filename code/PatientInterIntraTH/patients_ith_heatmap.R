
# load data
normalised.vsd <- readRDS(file='/data/losic_vst_norm_geneExp.rds')

normalised.vsd <- normalised.vsd[, c("H1.a", "H1.b", 
"H2.a1", "H2.b", "H2.c", "H2.d", "H2.e", #  "H2.a2"
"H3.a", "H3.b", 
"H4.a", "H4.b", "H4.c", "H4.d", "H4.e", 
"H6.a", "H6.b", 
"H7.a", "H7.b", "H7.c", "H7.d", "H7.e", 
"H8.a", "H8.b", "H8.c", 
"H9.a", "H9.b", "H9.c", "H9.d", "H9.e", "H9.f", 
"H10.a", "H10.b", "H10.c", "H10.d", "H10.e", 
"H11.a", "H11.b", 
"H12.a", "H12.b", "H12.c", "H12.d", "H12.e")]
colnames(normalised.vsd)[3] <- "H2.a"


# top 1000 variably expressed genes
exp.sd <- apply(normalised.vsd, 1, sd)
top.sd.genes <- names(exp.sd)[order(exp.sd, decreasing=TRUE)[1:1000]]
normalised.vsd <- normalised.vsd[top.sd.genes, ]


# clinical data
clinical.data <- data.frame(PatientID=paste0('H', c(1:4, 6:12)), 
 Etiology=c(rep('HBV', 4), rep('Cryptogenic', 2), 'HCV', 'Cryptogenic', 'HBV', 'HCV', 'HCV'))

# plot
IthHeatmap <- function(normalised.exp, patient.anno, file.name){
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(dendsort)
  library(ComplexHeatmap)
  library(viridis)
  library(circlize)
  

  # input: mat1, homoscedastic and library size-normalized count values
  mat1 <- normalised.exp %>% t() %>% as.data.frame() 
  
  # calculate Z-score
  tmp <- sapply(X=mat1, FUN=scale)
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- rownames(mat1)
  mat1 <- tmp
  
  # make mat2
  mat2 <- mat1 %>% rownames_to_column("SampleID")
  # add PatientID and RegionID cols
  mat2$PatientID <- sapply(X=mat2$SampleID, FUN=function(x) {unlist(strsplit(x, split="\\."))[1]})
  mat2$RegionID <- sapply(X=mat2$SampleID, FUN=function(x) {unlist(strsplit(x, split="\\."))[2]})
  # join clinical data
  mat2 <- dplyr::left_join(x = mat2, y=patient.anno, by="PatientID")

  mat2 <- dplyr::select(mat2, SampleID, PatientID, Etiology)
  # spread
  mat2 <- spread(mat2, SampleID, Etiology) 
  # transpose
  tmp <- mat2[,-1] %>% t() %>% as.data.frame()
  colnames(tmp) <- mat2$PatientID
  mat2 <- tmp
  
  # ComplexHeatmap functions
  mat2 <- mat2[rownames(mat1), ]
  
  # ht1: expression heatmap
  # reorder cols and rows
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  mat_cluster_cols <- hclust(dist(t(mat1)))
  mat_cluster_rows <- sort_hclust(hclust(dist(mat1)))
  # make heatmap-1
  ht1 = Heatmap(
    matrix=as.matrix(mat1), 
    name = "ht1", 
    # col = colorRamp2(seq(-2.7,4,1.5), viridis(5)),
	col = colorRamp2(seq(-3,3), viridis(7, option = "magma")),
    column_dend_height = unit(15, "mm"),
    row_dend_width = unit(15, "mm"),
    row_dend_reorder = rev(mat_cluster_rows$order), 
    column_dend_reorder = rev(mat_cluster_cols$order),
    show_column_names = FALSE,
    show_row_names = TRUE,
    width = 1,
    heatmap_legend_param = list(title = NULL, color_bar = "continuous")
  )
  

  # ht2: tumour region heatmap
  # make heatmap-2
  
  # specify colours
  colors = structure(c("green4", "darkorchid3", "darkorange"), names = c("HBV", "HCV", "Cryptogenic"))
 
  ht2 = Heatmap(
    matrix=as.matrix(mat2),
    col = colors,
    name = "ht2",
    na_col = "white",
    rect_gp = gpar(col="white"),
    show_row_names = TRUE, 
    show_column_names = TRUE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = NULL, color_bar = "discrete"),
    column_names_gp = gpar(fontsize = 4),
    width = 1
  )
  
  pdf(paste0(file.name, '.pdf'))
  # plot ComplexHeatmaps 
  print(ht1+ht2)

  #ht1: add border
  decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
  #ht2: add border, and lines
  decorate_heatmap_body("ht2", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
  dev.off()
}

setwd('/result/Section1')
IthHeatmap(normalised.vsd, clinical.data, 'losic_patient_ith_heatmap_1000sd')


##################################### Array
# Gene expression from Shi et al.
shi.gene.exp <- readRDS(file='/data/shi.lihc.gene.exp.rds')
shi.gene.exp <- shi.gene.exp[, -c(1:3)]

# top 1000 variably expressed genes
exp.sd <- apply(shi.gene.exp, 1, sd)
top.sd.genes <- names(exp.sd)[order(exp.sd, decreasing=TRUE)[1:1000]]
shi.gene.exp <- shi.gene.exp[top.sd.genes, ]


# clinical data
clinical.data <- data.frame(PatientID=paste0('H', 1:5), Etiology=c(rep('HBV', 5)))
IthHeatmap(shi.gene.exp, clinical.data, 'shi_patient_ith_heatmap_1000sd')




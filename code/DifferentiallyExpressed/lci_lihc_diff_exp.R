
# load data
cli.data <- readRDS(file='/data/lci_xinweiwang_cli_data.rds')
load(file='/data/lci_xinweiwang_expr_data.RData')
candi.prog.gene <- readRDS(file='/data/candi_prog_sig_gene.rds')


# data prepare
cli.data$ID <- gsub(' ', '', cli.data$ID)

tumor.id <- cli.data$ID[cli.data$Tissue.Type == 'Tumor']
normal.id <- cli.data$ID[cli.data$Tissue.Type == 'Non-Tumor']
paired.id <- intersect(tumor.id, normal.id)

tumor.sample <- subset(cli.data, Tissue.Type == 'Tumor' & ID %in% paired.id)$Affy_GSM
normal.sample <- subset(cli.data, Tissue.Type == 'Non-Tumor' & ID %in% paired.id)$Affy_GSM

eset <- gene.max.exp.profile[, c(tumor.sample, normal.sample)]

# Calculate differential expression genes
library(limma)
group <- rep(c('Tumor', 'Normal'), each=length(paired.id))
groupdesign <- model.matrix(~0+group)

rownames(groupdesign) <- c(tumor.sample, normal.sample)
colnames(groupdesign) <- c('Tumor', 'Normal')


fit <- lmFit(eset, groupdesign)
cont.matrix <- makeContrasts(TvsN=Normal-Tumor, levels = groupdesign)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


top.table <- topTable(fit2, sort.by = "P", n = Inf)
# top.table[candi.prog.gene, ]



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


setwd('/result/Section2/')
sam.anno <- data.frame(Type=c(rep('Tumor', length(tumor.sample)), rep('Normal', length(normal.sample))))
rownames(sam.anno) <- c(tumor.sample, normal.sample)

signature <- intersect(candi.prog.gene, rownames(eset))
ExpClusteringHeatmap(eset, signature, sam.anno, 'lci_candi_gene_exp_heatmap')



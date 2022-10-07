
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


# top 10000 variably expressed genes
exp.sd <- apply(normalised.vsd, 1, sd)
top.sd.genes <- names(exp.sd)[order(exp.sd, decreasing=TRUE)[1:10000]]
normalised.vsd <- normalised.vsd[top.sd.genes, ]


# Principal Components Analysis
pca.res <- prcomp(t(normalised.vsd), scale = TRUE)
pca.pc <- as.data.frame(pca.res$x)


# sample information
setwd('/data/OriginalData/Losic_Nat_Commun_2020_LiverMultiregionExp')
sam.info <- read.csv(file='sample_annotation.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
rownames(sam.info) <- sam.info$Sample

pca.pc <- merge(pca.pc, sam.info[, c('Patient', 'Etiology')], by='row.names')
pca.pc$Patient <- factor(pca.pc$Patient, levels=c(1:4, 5:12))

###################plot
library('ggpubr')
scatter.plot <- ggscatter(pca.pc, 'PC1', 'PC2', shape = 'Etiology', color = 'Patient', label = 'Row.names')

ggsave(scatter.plot, file = '/result/Section1/losic_pca10000.pdf')


##################################### Array
# Gene expression from Shi et al.
shi.gene.exp <- readRDS(file='/data/shi.lihc.gene.exp.rds')
shi.gene.exp <- shi.gene.exp[, -c(1:3)]

# top 10000 variably expressed genes
exp.sd <- apply(shi.gene.exp, 1, sd)
top.sd.genes <- names(exp.sd)[order(exp.sd, decreasing=TRUE)[1:length(exp.sd)]]
shi.gene.exp <- shi.gene.exp[top.sd.genes, ]


# Principal Components Analysis
pca.res <- prcomp(t(shi.gene.exp), scale = TRUE)
pca.pc <- as.data.frame(pca.res$x)


# clinical data
sam.info <- data.frame(Patient=rep(paste0('H', 1:5), each=5), Etiology=c(rep('HBV', 25)))
rownames(sam.info) <- paste0(rep(c('H1.', 'H2.', 'H3.', 'H4.', 'H5.'), each=5), 1:5)

pca.pc <- merge(pca.pc, sam.info[, c('Patient', 'Etiology')], by='row.names')
pca.pc$Patient <- factor(pca.pc$Patient, levels=paste0('H', 1:5))

###################plot
library('ggpubr')
scatter.plot <- ggscatter(pca.pc, 'PC1', 'PC2', shape = 'Etiology', color = 'Patient', label = 'Row.names')

ggsave(scatter.plot, file = '/result/Section1/shi_pca.pdf')



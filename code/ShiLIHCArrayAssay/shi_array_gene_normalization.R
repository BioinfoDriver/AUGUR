
# LIHC array dataset(GSE92528)

# Normalization
options(stringsAsFactors=F)
source('/code/Rscript/AgilentSCDataNormalization.R')

dir.path <- '/data/OriginalData/Shi_Oncotarget_2017_LiverMultiregionExp'
untar(tarfile = file.path(dir.path, 'GSE92528_RAW.tar'), exdir = file.path(dir.path, 'GSE92528'))


# #  Normalization of Agilent Single-Channel Data
file.path <- file.path(dir.path, 'GSE92528')
file.name = list.files(path = file.path, pattern = '.txt.gz')
gpl.num <- 'GPL6480'
gpl.anno.lab <- file.path(dir.path, 'Platforms')

# # keep probes that are above background on at least k arrays
# # # where k is the smallest number of replicates assigned to any of the experimental combinations
k.index <- 5
norm.exp.data <- AgilentSCDataNormalization(file.path, file.name, gpl.num, gpl.anno.lab, k.index)


# colnames
colnames(norm.exp.data)[-c(1:3)] <- do.call(rbind, 
 strsplit(colnames(norm.exp.data)[-c(1:3)], split = '_'))[, 1]
 
colnames(norm.exp.data)[-c(1:3)] <- paste0(rep(c('H1.', 'H2.', 'H3.', 'H4.', 'H5.'), each=5), 1:5)

# Protein coding gene
gene.info <- readRDS(file='/data/gene_info.rds')
norm.exp.data <- subset(norm.exp.data, GENE %in% gene.info$GeneID) # 15112

# save
saveRDS(norm.exp.data, file='/data/shi.lihc.gene.exp.rds')



# Accessing Raw Data from GEO
DownloadDataFromGEO <- function(file.path, acc.num){
 # GEO: A GEO accession number such as GSM1137.
 # baseDir: The base directory for the downloads.
 
 library('GEOquery')
 getGEOSuppFiles(GEO = acc.num, makeDirectory = TRUE, baseDir = file.path,
  fetch_files = TRUE, filter_regex = NULL)
 
 # Extract files from the contents of a tar archive
 untar(tarfile = paste0(c(file.path, acc.num, paste0(acc.num, '_RAW.tar')), collapse = '/'),
  exdir = paste0(c(file.path, acc.num), collapse = '/'))
}

# Agilent Single-Channel Data Normalization
AgilentSCDataNormalization <- function(file.path, file.name, gpl.num, gpl.anno.lab, k.index){
 # file.name: character vector giving the names of the files containing image analysis output
 # file.path: character string giving the directory containing the files.
 
 # k.index: keep probes that are above background on at least k arrays
 # where k is the smallest number of replicates assigned to any of the experimental combinations
 
 # Reading Single-Channel Agilent Intensity Data
 library('limma')
 exp.data <- read.maimages(files = file.name, source = 'agilent', 
  path = file.path, green.only=TRUE, other.columns = 'gIsWellAboveBG')

 # Background correction and normalize
 exp.data.bc <- backgroundCorrect(exp.data, method = 'normexp')
 exp.data.bcn <- normalizeBetweenArrays(exp.data.bc, method = 'quantile')
 
 # Gene annotation
 library('GEOquery')
 gpl.prob.anno <- getGEO(GEO = gpl.num, destdir = gpl.anno.lab)

 prob.anno <- gpl.prob.anno@dataTable@table
 prob.anno <- prob.anno[, c('ID', 'GENE', 'GENE_SYMBOL')]
 
 exp.data.bcn$genes <- cbind(exp.data.bcn$genes, 
  prob.anno[match(exp.data.bcn$genes$ProbeName, prob.anno$ID), c('GENE', 'GENE_SYMBOL')])
  
 # Gene filtering
 # Filter out control probes, those with no symbol, and those that don’t appear to be expressed
 Control <- exp.data.bcn$genes$ControlType %in% c(-1, 1)
 NoSymbol <-  exp.data.bcn$genes$GENE_SYMBOL == ''
 IsExpr <- rowSums(exp.data.bcn$other$gIsWellAboveBG > 0) >= k.index
 
 exp.data.bcn.filt <- exp.data.bcn[!Control & !NoSymbol & IsExpr, ]
 
 # For replicate probes, replace values with the mean
 exp.data.bcn.filt.mean <- avereps(exp.data.bcn.filt,
  ID = exp.data.bcn.filt$genes$ProbeName)
 
 # Select the probe with the highest normalized intensity averaged over all samples
 norm.probe.exp <- cbind(exp.data.bcn.filt.mean$genes[, c('GENE', 'GENE_SYMBOL', 'ProbeName')], 
  exp.data.bcn.filt.mean$E)
 
 norm.probe.exp <- split(x = norm.probe.exp, f = factor(norm.probe.exp$GENE))
 gene.exp.profile <- lapply(norm.probe.exp, function(exp.data){
  index <- which.max(rowSums(exp.data[, -c(1:3)]))
  exp.data <- exp.data[index, , FALSE]
 })
 
 gene.exp.profile <- do.call(rbind, gene.exp.profile)
 return(gene.exp.profile)
}


# demo

# Download the raw data with the specified GSM number from GEO
# dir.path <- '/pub5/xiaoyun/Jobs/J22/GEOExpressionData'
# gse.id <- 'GSE71187'
# DownloadDataFromGEO(file.path = dir.path, acc.num = gse.id)


# #  Normalization of Agilent Single-Channel Data
# file.path <- '/pub5/xiaoyun/Jobs/J22/GEOExpressionData/GSE71187'
# file.name = list.files(path = file.path, pattern = '.txt.gz')
# gpl.num <- 'GPL6480'
# gpl.anno.lab <- '/pub5/xiaoyun/Jobs/J22/GEOExpressionData/Platforms'

# # keep probes that are above background on at least k arrays
# # # where k is the smallest number of replicates assigned to any of the experimental combinations
# k.index <- 6
# norm.exp.data <- AgilentSCDataNormalization(file.path, file.name, gpl.num, gpl.anno.lab, k.index)

# # colnames(norm.exp.data)[-c(1:3)] <- do.call(rbind, 
 # # strsplit(colnames(norm.exp.data)[-c(1:3)], split = '_'))[, 1]

# save(norm.exp.data, file = 
 # '/pub5/xiaoyun/Jobs/J22/GEOExpressionData/NormalizedArrayExp/GSE71187.RData')



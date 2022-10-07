
# Accessing Clinical Data from GEO
DownloadCliDataFromGEO <- function(file.path, acc.num){
 # GEO: A GEO accession number such as GSM1137.
 # baseDir: The base directory for the downloads.
 
 library('GEOquery')
 gse <- getGEO(GEO = acc.num, destdir = file.path, GSEMatrix = TRUE)
 cli.data <- pData(phenoData(gse[[1]]))
 
 return(cli.data)
}


gse.id <- 'GSE54236'
dir.path <- 'data/OriginalData/GSE54236'
cli.data <- DownloadCliDataFromGEO(file.path = dir.path, acc.num = gse.id)

curated.cli.data <- cli.data[, c("title", "geo_accession", "gender:ch1", "doubling time (days):ch1", 
 "survival time(months):ch1", "tissue type:ch1")]
colnames(curated.cli.data) <- c("title", "geo_accession", "gender", "doubling_time", "survival_time", "tissue_type")

curated.cli.data$doubling_time <- as.numeric(curated.cli.data$doubling_time)
curated.cli.data$survival_time <- as.numeric(curated.cli.data$survival_time)


# Agilent Single-Channel Data Normalization
AgilentSCDataNormalization <- function(file.path, file.name, gpl.num, gpl.anno.lab, k.index, feac.meth){
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
 
 if(feac.meth=='max'){
  gene.exp.profile <- lapply(norm.probe.exp, function(exp.data){
   index <- which.max(rowSums(exp.data[, -c(1:3)]))
   exp.data <- exp.data[index, , FALSE]
  })
 }else{
  gene.exp.profile <- lapply(norm.probe.exp, function(exp.data){
   exp.data <- colMeans(exp.data[, -c(1:3)])
  }) 
 }

 gene.exp.profile <- do.call(rbind, gene.exp.profile)
 return(gene.exp.profile)
}



# #  Normalization of Agilent Single-Channel Data
file.path <- 'data/OriginalData/GSE54236/GSE54236_RAW'
file.name = list.files(path = file.path, pattern = '.txt.gz')
gpl.num <- 'GPL6480'
gpl.anno.lab <- 'data/OriginalData/GSE54236/Platforms'

# # keep probes that are above background on at least k arrays
# # # where k is the smallest number of replicates assigned to any of the experimental combinations
k.index <- 5
norm.exp.max.data <- AgilentSCDataNormalization(file.path, file.name, gpl.num, gpl.anno.lab, k.index, 'max')

colnames(norm.exp.max.data)[-c(1:3)] <- sapply(strsplit(colnames(norm.exp.max.data)[-c(1:3)], split = '_'), function(x) x[1])


save(norm.exp.max.data, file='/data/GSE54236_expr_data.RData')
saveRDS(curated.cli.data, file='/data/GSE54236_cli_data.rds')

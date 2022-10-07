
# Affymetrix microarray data normalization
# https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor
# http://www.bioconductor.org/packages/release/data/annotation/
# http://barc.wi.mit.edu/education/bioinfo2007/arrays/
# http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays
# Effects of filtering by Present call on analysis of microarray experiments

# Array normalization
AffymetrixArrayNormalization <- function(cel.path, file.names=NULL, probe.filter=FALSE, exp.threshold=5){
 library('affy')
 
 # Open CEL files from 3' Affymetrix Arrays (older ones) using affy
 if(is.null(file.names)){
  affy.data <- ReadAffy(celfile.path=cel.path)
 
 }else{
  affy.data <- ReadAffy(filenames=file.names, celfile.path=cel.path)
 
 }
 
 print(cdfName(affy.data))
 # Normalization using RMA
 eset.rma = rma(affy.data)
 expr.logs = exprs(eset.rma)


 # Filtering
 if(probe.filter){
  pa.calls  = mas5calls(affy.data)
  pa.calls = exprs(pa.calls)
  
  p.calls.n <- rowSums(pa.calls == 'P')
  expr.logs <- expr.logs[p.calls.n > exp.threshold, ]
 }

 return(expr.logs)
}


setwd('/data/OriginalData/GSE14520')
cli.data <- read.csv(file='GSE14520_Extra_Supplement.txt', header=T, sep='\t',stringsAsFactors=F)
GPL3921.anno <- GEOquery::getGEO(GEO='GPL3921', getGPL=FALSE, parseCharacteristics=FALSE)

GPL3921.samples <- intersect(GPL3921.anno@header$sample_id, cli.data$Affy_GSM)
cli.data <- subset(cli.data, Affy_GSM %in% GPL3921.samples)

cel.path <- '/data/OriginalData/GSE14520/GSE14520_RAW'
file.names <- paste0(GPL3921.samples, '.CEL.gz')
# at least 20% samples
norm.expr <- AffymetrixArrayNormalization(cel.path, file.names, probe.filter=TRUE, exp.threshold=89)

colnames(norm.expr) <- gsub('.CEL.gz', '', colnames(norm.expr))



# Probe annotation
probe.anno <- GPL3921.anno@dataTable@table
probe.anno <- subset(probe.anno, ENTREZ_GENE_ID != '')
probe.anno <- probe.anno[lengths(strsplit(probe.anno$ENTREZ_GENE_ID, split='///'))==1, ]
probe.anno <- probe.anno[, c('ID', 'ENTREZ_GENE_ID', 'Gene Symbol')]
rownames(probe.anno) <- probe.anno$ID

# Select the probe with the highest normalized intensity averaged over all samples
probe.norm.expr <- merge(probe.anno, norm.expr, by='row.names')
probe.norm.expr <- tibble::column_to_rownames(probe.norm.expr, var = 'Row.names')

probe.norm.expr <- split(x = probe.norm.expr, f = factor(probe.norm.expr$ENTREZ_GENE_ID))

# probe with maximal expression
gene.max.exp.profile <- lapply(probe.norm.expr, function(exp.data){
 index <- which.max(rowSums(exp.data[, -c(1:3)]))
 exp.data <- exp.data[index, , FALSE]
})

gene.max.exp.profile <- do.call(rbind, gene.max.exp.profile)


# save
save(gene.max.exp.profile, file='/data/lci_xinweiwang_expr_data.RData')
saveRDS(cli.data, file='/data/lci_xinweiwang_cli_data.rds')


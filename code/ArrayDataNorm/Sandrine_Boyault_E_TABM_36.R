
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


cel.path <- '/data/OriginalData/E-TABM-36/E-TABM-36.raw'
norm.expr <- AffymetrixArrayNormalization(cel.path, probe.filter=TRUE, exp.threshold=13)


# clinical data
cli.data <- read.csv(file='/data/OriginalData/E-TABM-36/E-TABM-36.sdrf.txt', 
 header=TRUE, sep='\t', stringsAsFactors=FALSE)

cli.data <- cli.data[, c(1, 45, 50, 51, 53:60, 62, 63)]
colnames(cli.data) <- c("PatientID", "Array", "Sex", "Age", "HBV_status", "Genotype", "DiseaseState",                 
"DiseaseStage", "FALIndex", "TP53_status", "os_status", "os_time", "pfs_status", "pfs_time")
rownames(cli.data) <- cli.data$PatientID

colnames(norm.expr) <- cli.data$PatientID[match(colnames(norm.expr), cli.data$Array)]

# Gene annotation
gpl.prob.anno <- getGEO(GEO = 'GPL96', destdir = '/data/OriginalData/E-TABM-36/Platforms')
probe.anno <- gpl.prob.anno@dataTable@table
probe.anno <- probe.anno[, c("ID", "Gene symbol", "Gene ID")]
colnames(probe.anno) <- c("ID", "GeneSymbol", "GeneID")


probe.anno <- subset(probe.anno, GeneID != '')
probe.anno <- probe.anno[lengths(strsplit(probe.anno$GeneID, split='///'))==1, ]
rownames(probe.anno) <- probe.anno$ID

# Select the probe with the highest normalized intensity averaged over all samples
probe.norm.expr <- merge(probe.anno, norm.expr, by='row.names')
probe.norm.expr <- tibble::column_to_rownames(probe.norm.expr, var = 'Row.names')

probe.norm.expr <- split(x = probe.norm.expr, f = factor(probe.norm.expr$GeneID))

# probe with maximal expression
gene.max.exp.profile <- lapply(probe.norm.expr, function(exp.data){
 index <- which.max(rowSums(exp.data[, -c(1:3)]))
 exp.data <- exp.data[index, , FALSE]
})

gene.max.exp.profile <- do.call(rbind, gene.max.exp.profile)


# save
save(gene.max.exp.profile, file='/data/E_TABM_36_expr_data.RData')
saveRDS(cli.data, file='/data/E_TABM_36_cli_data.rds')



setwd('/data/OriginalData/RefseqHsapiens')
gene.info <- read.csv(file='Homo_sapiens.gene_info', header = TRUE, sep = '\t', stringsAsFactors=FALSE)

# protein coding gene
gene.info <- subset(gene.info, type_of_gene == 'protein-coding')

gene.info$EnsemblID <- sapply(strsplit(gene.info$dbXrefs, split=':'), function(x)
	ifelse(length(grep('ENSG', x)) == 0, NA, x[grep('ENSG', x)]))

gene.info$EnsemblID <- sapply(strsplit(gene.info$EnsemblID, split='\\|'), function(x) x[1])
	

saveRDS(gene.info, file='/data/gene_info.rds')



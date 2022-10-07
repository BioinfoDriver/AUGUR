
#################GSE4024
library('GEOquery')
gse4024 <- getGEO(GEO = 'GSE4024', destdir = '/data/OriginalData/GSE4024and1898/GSE4024', GSEMatrix = TRUE)
exp.gse4024 <- exprs(gse4024[[1]])
cli.gse4024 <- pData(phenoData(gse4024[[1]]))
probe.anno.gse4024 <- gse4024[[1]]@featureData@data

# Clinical data
curated.cli.gse4024 <- cli.gse4024[, c("Individual:ch2", "geo_accession", "label_ch2", "Age:ch2", "Gender:ch2")]
colnames(curated.cli.gse4024) <- c("patient_id", "geo_accession", "label", "Age", "Gender")
pro.gse4024 <- read.csv(file='/data/OriginalData/GSE4024and1898/GSE4024/GSE4024.info.tsv', header=T, stringsAsFactors=F, sep='\t')
curated.cli.gse4024 <- merge(curated.cli.gse4024,  pro.gse4024, by.x='geo_accession', by.y='Array')


exp.gse4024 <- exp.gse4024[, subset(curated.cli.gse4024, label == 'cy5')$geo_accession] 


#################GSE1898
gse1898 <- getGEO(GEO = 'GSE1898', destdir = '/data/OriginalData/GSE4024and1898/GSE1898', GSEMatrix = TRUE)
cli.gse1898 <- pData(phenoData(gse1898[[1]]))

# Clinical data
curated.cli.gse1898 <- cli.gse1898[, c("title", "geo_accession")]
curated.cli.gse1898$patient_id  <- do.call(rbind, strsplit(curated.cli.gse1898$title, split='_'))[, 1]
curated.cli.gse1898$label  <- do.call(rbind, strsplit(curated.cli.gse1898$title, split='_'))[, 2]
pro.gse1898 <- read.csv(file='/data/OriginalData/GSE4024and1898/GSE1898/GSE1898.info.tsv', header=T, stringsAsFactors=F, sep='\t')
curated.cli.gse1898 <- merge(curated.cli.gse1898,  pro.gse1898, by.x='geo_accession', by.y='Array')


exp.gse1898 <- read.csv(file='/data/OriginalData/GSE4024and1898/GSE1898/GSE1898_series_matrix.txt', 
 sep='\t', comment.char = "!", header=TRUE, stringsAsFactors=FALSE)
exp.gse1898 <- tibble::column_to_rownames(exp.gse1898, var = "ID_REF")

exp.gse1898 <- exp.gse1898[, subset(curated.cli.gse1898, label == 'Cy5')$geo_accession] 

#################
# subset(probe.anno.gse4024, !(ID %in% rownames(exp.gse1898)))
       # ID  GENE UNIGENE                                                                    DESCRIPTION
# 1173521_1  FGL1  491143                          fibrinogen-like 1 (FGL1), transcript variant 1, mRNA.
# 1173531_1 CREB5  437075 cAMP responsive element binding protein 5 (CREB5), transcript variant 4, mRNA.
# 1173626_1   AFM  168718                                                            afamin (AFM), mRNA.
# 1173809_1  ALAD    1227       aminolevulinate, delta-, dehydratase (ALAD), transcript variant 2, mRNA.
# 1184076_1 TOP3A  592115                                   topoisomerase (DNA) III alpha (TOP3A), mRNA.
# 1184172_1  GYPE  632594                              glycophorin E (GYPE), transcript variant 2, mRNA.
# 1702056_1  ACTB  520640                                                                    Actin, beta
# 1702057_1  ACTB  520640                                                                    Actin, beta
# 1702058_1  ACTB  520640                                                                    Actin, beta
# 1702059_1  ACTB  520640                                                                    Actin, beta

#################
expr.data <- cbind(exp.gse1898, exp.gse4024[rownames(exp.gse1898), ])


# Batch rectification
library(sva)
batch.expr.data <- ComBat(dat = expr.data, batch = rep(c(1, 2), times=c(ncol(exp.gse1898), ncol(exp.gse4024))))


# Protein coding gene
{
probe.anno <- probe.anno.gse4024[, c('ID', 'GENE')] # 19746
# length(unique(probe.anno$GENE)), 16721
probe.anno <- subset(probe.anno, GENE!='' )

gene.info <- readRDS(file='/data/gene_info.rds')


# match symbol
symbol.gene.anno <- subset(gene.info, Symbol %in% probe.anno$GENE) # 11321
symbol.gene.anno$ExpGene <- symbol.gene.anno$Symbol

# match synonym
synonym.gene.anno <- subset(gene.info, Synonyms %in% probe.anno$GENE & !(Synonyms %in% Symbol))
synonym.gene.anno <- subset(synonym.gene.anno, !(Synonyms %in% Synonyms[duplicated(Synonyms)]))
synonym.gene.anno$ExpGene <- synonym.gene.anno$Synonyms # 410

# don't match
diff.gene <- setdiff(probe.anno$GENE, gene.info$Symbol)
diff.gene <- setdiff(diff.gene, gene.info$Synonyms) # 171


index <- which(is.na(match(gene.info$Symbol, probe.anno$GENE)) & is.na(match(gene.info$Synonyms, probe.anno$GENE)))
synonyms <- gene.info$Synonyms[index]
synonyms <- synonyms[lengths(strsplit(synonyms, split='\\|'))>1]
split.synonyms <- strsplit(synonyms, split='\\|')


diff.gene.anno <- lapply(diff.gene, function(gene.symbol){

	index.synonym <- which(!is.na(sapply(split.synonyms, function(x) match(gene.symbol, x))))
		
	if(length(index.synonym)>0){
		
		return(c(gene.symbol, synonyms[index.synonym[1]]))
	}else{
	
		return(c(gene.symbol, NA))
	}
	
})

diff.gene.anno <- do.call(rbind, diff.gene.anno)
colnames(diff.gene.anno) <- c('ExpGene', 'Synonyms')
diff.gene.anno <- as.data.frame(diff.gene.anno)

diff.gene.anno <- subset(diff.gene.anno, !is.na(Synonyms) & !duplicated(Synonyms))
diff.gene.anno <- merge(gene.info, diff.gene.anno, by = 'Synonyms') # 2418

exp.gene.anno <- rbind(symbol.gene.anno, synonym.gene.anno, diff.gene.anno[, colnames(symbol.gene.anno)])
exp.gene.anno <- subset(exp.gene.anno, !duplicated(Symbol)) # 14126

exp.gene.anno <- merge(exp.gene.anno, probe.anno, by.x='ExpGene', by.y='GENE') # 16843
rownames(exp.gene.anno) <- exp.gene.anno$ID
}

# Select the probe with the highest normalized intensity averaged over all samples
{
batch.expr.data <- merge(exp.gene.anno[, c('GeneID', 'Symbol', 'ID')], batch.expr.data, by='row.names')
batch.expr.data<- tibble::column_to_rownames(batch.expr.data, var = 'Row.names')

batch.expr.data <- split(x = batch.expr.data, f = factor(batch.expr.data$GeneID))
batch.expr.max.data <- lapply(batch.expr.data, function(exp.data){
 index <- which.max(rowSums(exp.data[, -c(1:3)], na.rm = TRUE))
 exp.data <- exp.data[index, , FALSE]
})

batch.expr.max.data <- do.call(rbind, batch.expr.max.data)

}

# save
save(batch.expr.max.data, file='/data/lec_snorri_s_thorgeirsson_expr_data.RData')

cli.data <- rbind(pro.gse1898, pro.gse4024)
cli.data$GSE <- rep(c('GSE1898', 'GSE4024'), times=c(nrow(pro.gse1898), nrow(pro.gse4024)))
saveRDS(cli.data, file='/data/lec_snorri_s_thorgeirsson_cli_data.rds')




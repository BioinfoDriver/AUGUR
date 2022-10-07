
# NCBI gene annotation
gene.info <- readRDS(file='/data/gene_info.rds')


# annotation of LOSIC's gene expression
load(file='/data/losic.refseq.gene.exp.RData')
losic.exp.gene <- do.call(rbind, strsplit(rownames(tran.exp$counts), split='_'))[, 2] # 19317

# match symbol
symbol.gene.anno <- subset(gene.info, Symbol %in% losic.exp.gene) # 19103
symbol.gene.anno$LosicExpGene <- symbol.gene.anno$Symbol

# match synonym
synonym.gene.anno <- subset(gene.info, Synonyms %in% losic.exp.gene & !(Synonyms %in% Symbol))
synonym.gene.anno <- subset(synonym.gene.anno, !(Synonyms %in% Synonyms[duplicated(Synonyms)]))
synonym.gene.anno$LosicExpGene <- synonym.gene.anno$Synonyms # 25

# don't match
diff.gene <- setdiff(losic.exp.gene, gene.info$Symbol)
diff.gene <- setdiff(diff.gene, gene.info$Synonyms) # 171


index <- which(is.na(match(gene.info$Symbol, losic.exp.gene)) & is.na(match(gene.info$Synonyms, losic.exp.gene)))
synonyms <- gene.info$Synonyms[index]
synonyms <- synonyms[lengths(strsplit(synonyms, split='\\|'))>1]
split.synonyms <- strsplit(synonyms, split='\\|')


diff.gene.anno <- lapply(diff.gene, function(gene.symbol){

	index.synonym <- which(!is.na(sapply(split.synonyms, function(x) match(gene.symbol, x))))
		
	if(length(index.synonym)>0){
		
		return(c(gene.symbol, synonyms[index.synonym]))
	}else{
	
		return(c(gene.symbol, NA))
	}
	
})

diff.gene.anno <- do.call(rbind, diff.gene.anno)
colnames(diff.gene.anno) <- c('LosicExpGene', 'Synonyms')
diff.gene.anno <- as.data.frame(diff.gene.anno)

diff.gene.anno <- subset(diff.gene.anno, !is.na(Synonyms) & !duplicated(Synonyms))
diff.gene.anno <- merge(gene.info, diff.gene.anno, by = 'Synonyms') # 43

losic.exp.gene.anno <- rbind(symbol.gene.anno, synonym.gene.anno, diff.gene.anno[, colnames(symbol.gene.anno)])
losic.exp.gene.anno <- subset(losic.exp.gene.anno, !duplicated(Symbol)) # 19170


# annotation of TCGA LIHC's gene expression
load(file='/data/tcga_lihc_rnaseq.RData')
tcga.exp.gene <- sapply(strsplit(rownames(tcga.lihc.count), split='\\|'), function(x) x[2]) # 20531
tcga.exp.gene.anno <- subset(gene.info, GeneID %in% tcga.exp.gene) # 18534
tcga.exp.gene.anno$TCGAExpGene <- tcga.exp.gene.anno$GeneID


# annotation of ICGC LIHC's gene expression
load(file='/data/icgc_lihc_filter_rnaseq.RData')
icgc.exp.gene <- rownames(icgc.lihc.filter.count) # 22913


# match symbol
symbol.gene.anno <- subset(gene.info, Symbol %in% icgc.exp.gene) # 16889
symbol.gene.anno$ICGCExpGene <- symbol.gene.anno$Symbol

# match synonym
synonym.gene.anno <- subset(gene.info, Synonyms %in% icgc.exp.gene & !(Synonyms %in% Symbol))
synonym.gene.anno <- subset(synonym.gene.anno, !(Synonyms %in% Synonyms[duplicated(Synonyms)]))
synonym.gene.anno$ICGCExpGene <- synonym.gene.anno$Synonyms # 408


# don't match
diff.gene <- setdiff(icgc.exp.gene, gene.info$Symbol)
diff.gene <- setdiff(diff.gene, gene.info$Synonyms) # 5614


index <- which(is.na(match(gene.info$Symbol, icgc.exp.gene)) & is.na(match(gene.info$Synonyms, icgc.exp.gene)))
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
colnames(diff.gene.anno) <- c('ICGCExpGene', 'Synonyms') 
diff.gene.anno <- as.data.frame(diff.gene.anno)

diff.gene.anno <- subset(diff.gene.anno, !is.na(Synonyms) & !duplicated(Synonyms))
diff.gene.anno <- merge(gene.info, diff.gene.anno, by = 'Synonyms') # 1345

icgc.exp.gene.anno <- rbind(symbol.gene.anno, synonym.gene.anno, diff.gene.anno[, colnames(symbol.gene.anno)])
icgc.exp.gene.anno <- subset(icgc.exp.gene.anno, !duplicated(Symbol)) # 18634

# combine
exp.gene.anno <- merge(gene.info, tcga.exp.gene.anno[, c('GeneID', 'TCGAExpGene')], by='GeneID', all.x=T)
exp.gene.anno <- merge(exp.gene.anno, icgc.exp.gene.anno[, c('GeneID', 'ICGCExpGene')], by='GeneID', all.x=T)
exp.gene.anno <- merge(exp.gene.anno, losic.exp.gene.anno[, c('GeneID', 'LosicExpGene')], by='GeneID', all.x=T)

exp.gene.anno <- exp.gene.anno[, c(2, 1, 3:20)]
save(exp.gene.anno, file='/data/exp_gene_anno.RData') # 18339

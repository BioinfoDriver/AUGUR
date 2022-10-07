
# load data
tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')

# Published signature
library('xlsx')
lihc.signature <- read.xlsx(file='/data/Signature.xlsx', sheetIndex =1)

library(org.Hs.eg.db)
lihc.signature$GeneID <- mapIds(org.Hs.eg.db, keys = lihc.signature$GeneName, keytype = "SYMBOL", column="ENTREZID")
lihc.signature <- subset(lihc.signature, !is.na(GeneID)) 
# Lim_AnnSurgOncol_2013(23800896) 少2个基因ILMN_1851092, ILMN_1868912

# saveRDS(lihc.signature, file='/data/curated_lihc_risk_signature.rds')

lihc.signature <- subset(lihc.signature, GeneID %in% rownames(tcga.lihc.vst.exp))
# Li_FrontCellDevBiol_2022(35047511) 少1个基因——PET117
# Zhang_JCellPhysiol_2019(30556603) 少1个基因——PGM5P3-AS1


# risk score
signature.risk.score <- sapply(unique(lihc.signature$PMID), function(ID){
	
	risk.coef <- subset(lihc.signature, PMID == ID)
	
	exp.dat <- as.matrix(tcga.lihc.vst.exp[risk.coef$GeneID, ])
	risk.score <- crossprod(exp.dat, as.matrix(risk.coef$Coefficient))[, 1]
	
	return(risk.score)
})

colnames(signature.risk.score) <- paste0('PMID:', unique(lihc.signature$PMID))

saveRDS(signature.risk.score, file='/data/curated_sig_tcga_risk_score.rds')



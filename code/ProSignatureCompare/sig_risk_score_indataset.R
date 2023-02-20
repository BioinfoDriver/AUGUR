
# gene signatures
lihc.signature <- readRDS(file='/data/curated_lihc_risk_signature.rds')
load(file='/data/lasso_cox_res.RData')

lihc.signature <- lihc.signature[, c('GeneID', 'Coefficient', 'PMID')]
lihc.signature <- subset(lihc.signature, PMID %in% c('32198063', '30885723', '35123387', '31335995', '31695766',
 '33033585', '34676211', '34277622', '33777771', '34975331', '25666192','22105560'))

lihc.signature <- rbind(data.frame(GeneID = names(active.k.vals),  
 Coefficient = active.k.vals, PMID = 'AUGUR'), lihc.signature)

nault.sig.class.label <- readRDS(file='/data/nault_sig_class_label.rds')

# function

ProSignatureScore <- function(signature, nault.sig=NULL, exp.data, cli.data){

	# intersect
	signature <- subset(signature, GeneID %in% rownames(exp.data))
	
	
	# risk score
	signature.risk.score <- sapply(unique(signature$PMID), function(ID){
		
		risk.coef <- subset(signature, PMID == ID)
		
		exp <- as.matrix(exp.data[risk.coef$GeneID, ])
		risk.score <- crossprod(exp, as.matrix(risk.coef$Coefficient))[, 1]
		
		if(ID == '33777771'){
			risk.score <- exp(risk.score)
		}
		
		return(risk.score)	
	})
	
	signature.risk.score <- as.data.frame(signature.risk.score)
	colnames(signature.risk.score) <- paste0('PMID_', unique(signature$PMID))
	
	if(!is.null(nault.sig)){
		
		nault.sig <- ifelse(nault.sig == 'Bad', 1, 0)
		signature.risk.score$PMID_23567350 <- nault.sig[rownames(signature.risk.score)]
	
	}
	
	cli.data <- cbind(cli.data, signature.risk.score[rownames(cli.data), ])
	
	return(cli.data)
}


# load data

# TCGA
tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')
tcga.lihc.cli.data <- readRDS(file='/data/tcga.lihc.cli.data.rds')
rownames(tcga.lihc.cli.data) <- paste0(rownames(tcga.lihc.cli.data), '-01')
tcga.lihc.cli.data <- tcga.lihc.cli.data[colnames(tcga.lihc.vst.exp), ]
tcga.lihc.cli.data <- tcga.lihc.cli.data[, c('os', 'os_time')]
colnames(tcga.lihc.cli.data) <- c('os', 'os.time')

nault.sig <- nault.sig.class.label$TCGA_LIHC
names(nault.sig) <- paste0(names(nault.sig), '-01')

sig.risk.score.in.tcga <- ProSignatureScore(lihc.signature, nault.sig, tcga.lihc.vst.exp, tcga.lihc.cli.data)


# ICGC
icgc.lihc.vst.exp <- readRDS(file='/data/icgc_vst_norm_geneExp.rds')
icgc.linc.cli.data <- readRDS(file='/data/icgc_linc_cli_data.rds')
icgc.linc.cli.data <- subset(icgc.linc.cli.data, !filter)
rownames(icgc.linc.cli.data) <- icgc.linc.cli.data$icgc_sample_id

icgc.linc.cli.data <- icgc.linc.cli.data[, c('donor_vital_status', 'donor_survival_time')]
colnames(icgc.linc.cli.data) <- c('os', 'os.time')
icgc.linc.cli.data$os <- ifelse(icgc.linc.cli.data$os == 'deceased', 1, 0)

sig.risk.score.in.icgc <- ProSignatureScore(lihc.signature, 
 nault.sig.class.label$ICGC_LIRI, icgc.lihc.vst.exp, icgc.linc.cli.data)


# CHCC-HBV
CHCC.HBV.norm.exp <- readRDS(file='/data/Gao_Cell_2019_expr_data.rds')
CHCC.HBV.cli.data <- readRDS(file='/data/Gao_Cell_2019_cli_data.rds')

rownames(CHCC.HBV.norm.exp) <- substr(rownames(CHCC.HBV.norm.exp), 1, 15)
library(org.Hs.eg.db)
rownames(CHCC.HBV.norm.exp) <- mapIds(org.Hs.eg.db, keys = rownames(CHCC.HBV.norm.exp), keytype = "ENSEMBL", column="ENTREZID")

CHCC.HBV.cli.data <- CHCC.HBV.cli.data[, c("Survial  (1, dead; 0, alive)", "Overall survial (month)")]
colnames(CHCC.HBV.cli.data) <- c('os', 'os.time')

sig.risk.score.in.chcc <- ProSignatureScore(lihc.signature, 
 nault.sig.class.label$CHCC_HBV, CHCC.HBV.norm.exp, CHCC.HBV.cli.data)


# Mongolian-HCC
Mongolian.HCC.norm.exp <- readRDS(file='/data/GSE144269_expr_data.rds')
Mongolian.HCC.cli.data <- readRDS(file='/data/GSE144269_cli_data.rds')
rownames(Mongolian.HCC.cli.data) <- Mongolian.HCC.cli.data$RNASeq_T

Mongolian.HCC.cli.data <- Mongolian.HCC.cli.data[, c('survival.status', 'survival.time')]
colnames(Mongolian.HCC.cli.data) <- c('os', 'os.time')

rownames(Mongolian.HCC.norm.exp) <- substr(rownames(Mongolian.HCC.norm.exp), 1, 15)
rownames(Mongolian.HCC.norm.exp) <- mapIds(org.Hs.eg.db, 
 keys = rownames(Mongolian.HCC.norm.exp), keytype = "ENSEMBL", column="ENTREZID")

sig.risk.score.in.mongolian <- ProSignatureScore(lihc.signature, nault.sig.class.label$Mongolian_HCC,
 Mongolian.HCC.norm.exp, Mongolian.HCC.cli.data)
 

# FULCI-HCC
lci.cli.data <- readRDS(file='/data/lci_xinweiwang_cli_data.rds')
load(file='/data/lci_xinweiwang_expr_data.RData')
lci.exp.data <- gene.max.exp.profile

lci.cli.data <- subset(lci.cli.data, Tissue.Type == 'Tumor' & !is.na(Survival.status))
rownames(lci.cli.data) <- lci.cli.data$Affy_GSM
lci.cli.data <- lci.cli.data[, c('Survival.status', 'Survival.months')]
colnames(lci.cli.data) <- c('os', 'os.time')

# transform labels from months to days
lci.cli.data$os.time <- lci.cli.data$os.time*30.4375
lci.exp.data <- lci.exp.data[, -c(1:3)]

sig.risk.score.in.fulci <- ProSignatureScore(lihc.signature, nault.sig.class.label$FULCI_HCC, lci.exp.data, lci.cli.data)
 

# NCI-HCC
lec.cli.data <- readRDS(file='/data/lec_snorri_s_thorgeirsson_cli_data.rds')
load('/data/lec_snorri_s_thorgeirsson_expr_data.RData')
lec.exp.data <- batch.expr.max.data

# transform labels from months to days
lec.cli.data <- subset(lec.cli.data, !is.na(OS_Status))
lec.cli.data$OS_Time <- lec.cli.data$OS_Time*30.4375
rownames(lec.cli.data) <- lec.cli.data$Array
lec.cli.data <- lec.cli.data[, c('OS_Status', 'OS_Time')]
colnames(lec.cli.data) <- c('os', 'os.time')
 
lec.exp.data <- lec.exp.data[, -c(1:3)]
lec.exp.data <- lec.exp.data[, intersect(rownames(lec.cli.data), colnames(lec.exp.data))]
lec.cli.data <- lec.cli.data[colnames(lec.exp.data), ]

sig.risk.score.in.nci <- ProSignatureScore(lihc.signature, nault.sig.class.label$NCI_HCC, lec.exp.data, lec.cli.data)
 

# INSERM-HCC
INSERM.HCC.cli.data <- readRDS(file='/data/E_TABM_36_cli_data.rds')
load(file='/data/E_TABM_36_expr_data.RData')
INSERM.HCC.exp.data <- gene.max.exp.profile

INSERM.HCC.cli.data <- subset(INSERM.HCC.cli.data, !is.na(os_status) & DiseaseState %in% c("HCC tumor"))
INSERM.HCC.cli.data$os_time <- INSERM.HCC.cli.data$os_time*30.4375
INSERM.HCC.cli.data <- INSERM.HCC.cli.data[, c('os_status', 'os_time')] 
colnames(INSERM.HCC.cli.data) <- c('os', 'os.time')
 
INSERM.HCC.exp.data <- INSERM.HCC.exp.data[, -c(1:3)]
INSERM.HCC.exp.data <- INSERM.HCC.exp.data[, rownames(INSERM.HCC.cli.data)]

sig.risk.score.in.inserm <- ProSignatureScore(lihc.signature, 
 nault.sig.class.label$INSERM_HCC, INSERM.HCC.exp.data, INSERM.HCC.cli.data)



# save data
sig.risk.score.indataset <- list('TCGA_LIHC'=sig.risk.score.in.tcga, 'ICGC_LIRI'=sig.risk.score.in.icgc, 
 'CHCC_HBV'=sig.risk.score.in.chcc, 'Mongolian_HCC'=sig.risk.score.in.mongolian, 
 'FULCI_HCC'=sig.risk.score.in.fulci, 'NCI_HCC'=sig.risk.score.in.nci, 'INSERM_HCC'=sig.risk.score.in.inserm)

saveRDS(sig.risk.score.indataset, file='/data/sig_risk_score_indataset.rds')


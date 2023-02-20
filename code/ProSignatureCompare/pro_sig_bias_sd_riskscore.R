
# load data
# Gene expression from Losic et al.
losic.gene.exp <- readRDS(file='/data/losic_vst_norm_geneExp.rds')


losic.gene.exp <- losic.gene.exp[, c(
"H1.a", "H1.b", 
"H2.a1", "H2.b", "H2.c", "H2.d", "H2.e", #  "H2.a2"
"H3.a", "H3.b", 
"H4.a", "H4.b", "H4.c", "H4.d", "H4.e", 
"H6.a", "H6.b", 
"H7.a", "H7.b", "H7.c", "H7.d", "H7.e", 
"H8.a", "H8.b", "H8.c", 
"H9.a", "H9.b", "H9.c", "H9.d", "H9.e", "H9.f", 
"H10.a", "H10.b", "H10.c", "H10.d", "H10.e", 
"H11.a", "H11.b", 
"H12.a", "H12.b", "H12.c", "H12.d", "H12.e")]
colnames(losic.gene.exp)[3] <- "H2.a"


# Gene expression from Shi et al.
shi.gene.exp <- readRDS(file='/data/shi.lihc.gene.exp.rds')
shi.gene.exp <- shi.gene.exp[, c(
"H1.1", "H1.2", "H1.3", "H1.4", "H1.5",
"H2.1", "H2.2", "H2.3", "H2.4", "H2.5", 
"H3.1", "H3.2", "H3.3", "H3.4", "H3.5", 
"H4.1", "H4.2", "H4.3", "H4.4", "H4.5", 
"H5.1", "H5.2", "H5.3", "H5.4", "H5.5")]


RiskScoreSD <- function(exp.data, risk.sig, sig.name){
    
  # input: signature 
  # input: expression matrix

  
  # select signature genes
  inter.sect <- intersect(rownames(exp.data), names(risk.sig))
  
  if(length(risk.sig) > length(inter.sect)){
	print(sig.name)
	print('The following genes are missed:')
	
	print(setdiff(names(risk.sig), inter.sect))
  }
  
  risk.sig <- risk.sig[inter.sect]
  exp.data <- exp.data[inter.sect, ]

 
  # calculate RiskScore
  RiskScore <- colSums(exp.data * risk.sig)
  
  if(sig.name == 'PMID33777771'){
	RiskScore <- exp(RiskScore)
  }
  
  RiskScore <- data.frame(SampleID=names(RiskScore), RiskScore=RiskScore)
  RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[1]})
  
  # SD
  risk.score.sd <- tapply(RiskScore$RiskScore, factor(RiskScore$PatientID), sd)
  return(risk.score.sd)
}

SigScoreSDPlot <- function(risk.sig.list, gene.exp.mat, file.name){
  
  risk.score.sds <- lapply(names(risk.sig.list), function(sig.name){
    risk.score.sd <- RiskScoreSD(gene.exp.mat, risk.sig.list[[sig.name]], sig.name)

    return(risk.score.sd)
  })

  sig.risk.score.sd <- data.frame(risk.score.sd=unlist(risk.score.sds), 
   sig.name = rep(names(risk.sig.list), times=sapply(risk.score.sds, length)))
  
  sig.risk.score.sd$sig.name <- factor(sig.risk.score.sd$sig.name, 
   levels=names(risk.sig.list)[order(sapply(risk.score.sds, median))])
  
  # print(rank(sapply(risk.score.sds, median)))
  
  # plot 
  library('ggpubr')
  score.sd.plot <- ggboxplot(sig.risk.score.sd, x='sig.name', y='risk.score.sd', fill='sig.name', 
   xlab='Signature', ylab='SD of risk score', legend = "none", add='point', x.text.angle=30)
  
  ggsave(score.sd.plot, file=paste(file.name, '.pdf'), width=7)
}


# Signature
# signature coefficient
load(file='/data/lasso_cox_res.RData')
sig.coef <- readRDS(file='/data/curated_lihc_risk_signature.rds')

risk.sig.list <- lapply(unique(sig.coef$PMID), function(pmid){
  
  risk.sig <- subset(sig.coef, PMID == pmid)$Coefficient
  names(risk.sig) <- subset(sig.coef, PMID == pmid)$GeneID
  
  return(risk.sig)
})
names(risk.sig.list) <- paste0('PMID', unique(sig.coef$PMID))
risk.sig.list <- c(list('CloSignature'=active.k.vals), risk.sig.list)

# RiskScoreSD(losic.gene.exp, risk.sig.list[['Clonal Signature']])

risk.sig <- c('CloSignature', 'PMID32198063', 'PMID30885723', 'PMID35123387', 'PMID31335995', 'PMID31695766',
 'PMID33033585', 'PMID34676211', 'PMID34277622', 'PMID33777771', 'PMID34975331', 'PMID25666192','PMID22105560')
risk.sig.list <- risk.sig.list[risk.sig] 

setwd('/result/Section5')
SigScoreSDPlot(risk.sig.list, losic.gene.exp, 'sig_score_sd_in_losic')
SigScoreSDPlot(risk.sig.list, shi.gene.exp, 'sig_score_sd_in_shi')


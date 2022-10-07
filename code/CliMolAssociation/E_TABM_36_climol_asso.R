
# load data
cli.data <- readRDS(file='/data/E_TABM_36_cli_data.rds')
load(file='/data/E_TABM_36_expr_data.RData')
load(file='/data/lasso_cox_res.RData')


norm.exp <- gene.max.exp.profile
hcc.cli.data <- cli.data
hcc.cli.data <- subset(hcc.cli.data, DiseaseState %in% c("HCC tumor"))
norm.exp <- norm.exp[, rownames(hcc.cli.data)]


# risk evaluation
RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
	
	# exp.dat <- as.matrix(exp.dat[intersect(rownames(exp.dat), gene.set), ])
	exp.dat <- as.matrix(exp.dat[gene.set, ])
	
	risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
	if(!is.null(cut.off)){
		risk.categ <- ifelse(risk.score >= cut.off, 'high risk', 'low risk')
	
	}else{
		risk.categ <- ifelse(risk.score >= median(risk.score), 'high risk', 'low risk')
		
	}
	return(data.frame(risk.score, risk.categ))
}


active.k.vals <- active.k.vals[intersect(names(active.k.vals), rownames(norm.exp))]# 11/12 
risk.score <- RiskEsti(exp.dat=norm.exp, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)
hcc.cli.data <- cbind(hcc.cli.data, risk.score[rownames(hcc.cli.data), ])



hcc.cli.data$Age <- ifelse(hcc.cli.data$Age < 60, '<60', '≥60')
hcc.cli.data$HBV_status <- ifelse(hcc.cli.data$HBV_status == 'HBV_titer_negative', 'No', 'Yes')
hcc.cli.data$Genotype[hcc.cli.data$Genotype == ''] <- NA
hcc.cli.data$Genotype <- factor(hcc.cli.data$Genotype, c('beta-catenine_not_mutated', 'beta-catenine_mutated'))
hcc.cli.data$DiseaseStage[hcc.cli.data$DiseaseStage == '  '] <- NA

hcc.cli.data$TP53_status <- factor(hcc.cli.data$TP53_status, levels=c('TP53_not_mutated', 'TP53_mutated'))
hcc.cli.data$risk.categ <- factor(hcc.cli.data$risk.categ, levels=c('low risk', 'high risk'))


# Fisher's Exact Test
cli.sig.char <- hcc.cli.data
p.values <- lapply(colnames(cli.sig.char)[c(4, 3, 5, 6, 8, 10)], function(char){
 print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
 print(table(cli.sig.char[, c(char, 'risk.categ')]))
 
 or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
 p.value <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
 return(c(or, p.value))
})

p.values <- as.data.frame(do.call(rbind, p.values))
colnames(p.values) <- c('OR', 'Pvalue')
p.values$Qvalue <- p.adjust(p.values$Pvalue, 'fdr')
rownames(p.values) <- colnames(cli.sig.char)[c(4, 3, 5, 6, 8, 10)]




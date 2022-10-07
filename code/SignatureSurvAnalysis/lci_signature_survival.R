
# load data
cli.data <- readRDS(file='/data/lci_xinweiwang_cli_data.rds')
load(file='/data/lci_xinweiwang_expr_data.RData')
load(file='/data/lasso_cox_res.RData')


lci.exp.data <- gene.max.exp.profile
lci.cli.data <- cli.data
# data prepare
lci.cli.data <- subset(lci.cli.data, Tissue.Type == 'Tumor')
lci.cli.data <- subset(lci.cli.data, !is.na(Survival.status))

# transform labels from months to days
lci.cli.data$Survival.months <- lci.cli.data$Survival.months*30.4375
lci.cli.data$Recurr.months <- lci.cli.data$Recurr.months*30.4375

lci.exp.data <- lci.exp.data[, lci.cli.data$Affy_GSM]


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

active.k.vals <- active.k.vals[intersect(names(active.k.vals), rownames(lci.exp.data))]# 12/12
lci.risk.score <- RiskEsti(exp.dat=lci.exp.data, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)
lci.cli.data <- cbind(lci.cli.data, lci.risk.score[lci.cli.data$Affy_GSM, ])



# survival plot 
source('/code/Rscript/survival_plot.R')
SurvivalPlot(survival.data=lci.cli.data[, c('LCS.ID', 'Survival.months', 'Survival.status')], 
 sample.class=lci.cli.data[, c('LCS.ID', 'risk.categ')], filename='lci_lihc_max_os.pdf', 
 out.file.path='/result/Section3/')

SurvivalPlot(survival.data=lci.cli.data[, c('LCS.ID', 'Recurr.months', 'Recurr.status')], 
 sample.class=lci.cli.data[, c('LCS.ID', 'risk.categ')], filename='lci_lihc_max_rfs.pdf', 
 out.file.path='/result/Section3/')


saveRDS(lci.cli.data, file='/data/lci_risk_score.rds')

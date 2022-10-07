
# load data
cli.data <- readRDS(file='/data/E_TABM_36_cli_data.rds')
load(file='/data/E_TABM_36_expr_data.RData')
load(file='/data/lasso_cox_res.RData')


norm.exp <- gene.max.exp.profile
hcc.cli.data <- cli.data

hcc.cli.data$os_time <- hcc.cli.data$os_time*30.4375
hcc.cli.data$pfs_time <- hcc.cli.data$pfs_time*30.4375
hcc.cli.data <- subset(hcc.cli.data, !is.na(os_status) & DiseaseState %in% c("HCC tumor"))
# hcc.cli.data <- subset(hcc.cli.data, DiseaseState %in% c("HCC tumor"))

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


# survival plot 
source('/code/Rscript/survival_plot.R')
SurvivalPlot(survival.data=hcc.cli.data[, c('PatientID', 'os_time', 'os_status')], 
 sample.class=hcc.cli.data[, c('PatientID', 'risk.categ')], filename='E_TABM_36_max_os.pdf', 
 out.file.path='/result/Section3/')

SurvivalPlot(survival.data=hcc.cli.data[, c('PatientID', 'pfs_time', 'pfs_status')], 
 sample.class=hcc.cli.data[, c('PatientID', 'risk.categ')], filename='E_TABM_36_max_pfs.pdf', 
 out.file.path='/result/Section3/')


saveRDS(hcc.cli.data, file='/data/E_TABM_36_risk_score.rds')

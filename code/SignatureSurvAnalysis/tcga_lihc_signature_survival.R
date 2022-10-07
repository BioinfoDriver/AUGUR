
# load data
load(file='/data/lasso_cox_res.RData')
tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')
tcga.lihc.cli.data <- readRDS(file='/data/tcga.lihc.cli.data.rds')


# data prepare
# active.k.vals <- active.k.vals.1se
colnames(tcga.lihc.vst.exp) <- substr(colnames(tcga.lihc.vst.exp), 1, 12)


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

tcga.risk.score <- RiskEsti(exp.dat=tcga.lihc.vst.exp, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)
tcga.lihc.cli.data <- merge(tcga.lihc.cli.data, tcga.risk.score, by='row.names')
colnames(tcga.lihc.cli.data)[1] <- 'patient_id'


# survival plot 
source('/code/Rscript/survival_plot.R')
SurvivalPlot(survival.data=tcga.lihc.cli.data[, c('patient_id', 'os_time', 'os')], 
 sample.class=tcga.lihc.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lihc_os.pdf', 
 out.file.path='/result/Section3/')

SurvivalPlot(survival.data=tcga.lihc.cli.data[, c('patient_id', 'dss_time', 'dss')], 
 sample.class=tcga.lihc.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lihc_dss.pdf', 
 out.file.path='/result/Section3/')

SurvivalPlot(survival.data=tcga.lihc.cli.data[, c('patient_id', 'pfi_time', 'pfi')], 
 sample.class=tcga.lihc.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lihc_pfi.pdf', 
 out.file.path='/result/Section3/')

SurvivalPlot(survival.data=tcga.lihc.cli.data[, c('patient_id', 'dfi_time', 'dfi')], 
 sample.class=tcga.lihc.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lihc_dfi.pdf', 
 out.file.path='/result/Section3/')


saveRDS(tcga.lihc.cli.data, file='/data/tcga_risk_score.rds')


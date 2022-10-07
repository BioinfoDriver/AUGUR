
# load data
load(file='/data/lasso_cox_res.RData')
icgc.linc.cli.data <- readRDS(file='/data/icgc_linc_cli_data.rds')
icgc.lihc.vst.exp <- readRDS(file='/data/icgc_vst_norm_geneExp.rds')


# data prepare
# active.k.vals <- active.k.vals.1se
icgc.linc.cli.data <- subset(icgc.linc.cli.data, !filter)
rownames(icgc.linc.cli.data) <- icgc.linc.cli.data$icgc_sample_id


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


# all patients
icgc.risk.score <- RiskEsti(exp.dat=icgc.lihc.vst.exp, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)
icgc.linc.cli.data <- merge(icgc.linc.cli.data, icgc.risk.score, by='row.names')



# survival plot 
source('/code/Rscript/survival_plot.R')
icgc.linc.cli.data$donor_vital_status <- ifelse(icgc.linc.cli.data$donor_vital_status == 'deceased', 1, 0)
SurvivalPlot(survival.data=icgc.linc.cli.data[, c('icgc_donor_id', 'donor_survival_time', 'donor_vital_status')], 
 sample.class=icgc.linc.cli.data[, c('icgc_donor_id', 'risk.categ')], filename='icgc_lihc_os.pdf', 
 out.file.path='/result/Section3/')

saveRDS(icgc.linc.cli.data, file='/data/icgc_risk_score.rds')


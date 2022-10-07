
# load data
cli.data <- readRDS(file='/data/lec_snorri_s_thorgeirsson_cli_data.rds')
load('/data/lec_snorri_s_thorgeirsson_expr_data.RData')
load(file='/data/lasso_cox_res.RData')

# data prepare
lec.exp.data <- batch.expr.max.data
lec.cli.data <- subset(cli.data, !is.na(OS_Status))
lec.exp.data <- lec.exp.data[, intersect(lec.cli.data$Array, colnames(lec.exp.data))]
lec.cli.data <- subset(lec.cli.data, Array %in% colnames(lec.exp.data))


# transform labels from months to days
lec.cli.data$OS_Time <- lec.cli.data$OS_Time*30.4375
rownames(lec.cli.data) <- lec.cli.data$Array

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

active.k.vals <- active.k.vals[intersect(names(active.k.vals), rownames(lec.exp.data))] #10/12
lec.risk.score <- RiskEsti(exp.dat=lec.exp.data, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)
lec.cli.data <- cbind(lec.cli.data, lec.risk.score[lec.cli.data$Array, ])




# survival plot 
source('/code/Rscript/survival_plot.R')
SurvivalPlot(survival.data=lec.cli.data[, c('Array', 'OS_Time', 'OS_Status')], 
 sample.class=lec.cli.data[, c('Array', 'risk.categ')], filename='lec_batch_max_os.pdf', out.file.path='/data/Section3/')

saveRDS(lec.cli.data, file='/data/lec_risk_score.rds')

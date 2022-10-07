
# load data
cli.data <- readRDS(file='/data/GSE54236_cli_data.rds')
load(file='/data/GSE54236_expr_data.RData')
load(file='/data/lasso_cox_res.RData')


GSE54236.exp.data <- norm.exp.max.data
GSE54236.cli.data <- cli.data

# data prepare
GSE54236.cli.data <- subset(GSE54236.cli.data, tissue_type == 'Biopsy of tumor tissue')
GSE54236.cli.data <- GSE54236.cli.data[-grep('Tumor_rep1', GSE54236.cli.data$title),]
# GSE54236.cli.data <- subset(GSE54236.cli.data, !(geo_accession %in% c('GSM1310584', 'GSM1310604', 'GSM1310610')))

# transform labels from months to days
GSE54236.cli.data$survival_time <- GSE54236.cli.data$survival_time*30.4375
GSE54236.cli.data$survival_status <- 1
GSE54236.cli.data$doubling_status <- 1

GSE54236.exp.data <- GSE54236.exp.data[, GSE54236.cli.data$geo_accession]


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

active.k.vals <- active.k.vals[intersect(names(active.k.vals), rownames(GSE54236.exp.data))]# 12/12, 12/12(filter)
GSE54236.risk.score <- RiskEsti(exp.dat=GSE54236.exp.data, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)
GSE54236.cli.data <- cbind(GSE54236.cli.data, GSE54236.risk.score[GSE54236.cli.data$geo_accession, ])


# survival plot 
source('/code/Rscript/survival_plot.R')
SurvivalPlot(survival.data=GSE54236.cli.data[, c('geo_accession', 'survival_time', 'survival_status')], 
 sample.class=GSE54236.cli.data[, c('geo_accession', 'risk.categ')], filename='GSE54236_expr_max_os.pdf', 
 out.file.path='/result/Section3/')

SurvivalPlot(survival.data=GSE54236.cli.data[, c('geo_accession', 'doubling_time', 'doubling_status')], 
 sample.class=GSE54236.cli.data[, c('geo_accession', 'risk.categ')], filename='GSE54236_expr_max_dbt.pdf', 
 out.file.path='/result/Section3/')
 

saveRDS(GSE54236.cli.data, file='/data/GSE54236_risk_score.rds')


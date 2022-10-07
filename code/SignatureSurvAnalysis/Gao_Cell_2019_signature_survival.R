
# load data
load(file='/data/lasso_cox_res.RData')
norm.exp <- readRDS(file='/data/Gao_Cell_2019_expr_data.rds')
cli.data <- readRDS(file='/data/Gao_Cell_2019_cli_data.rds')


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


active.k.vals <- c(-0.020624860,0.033253048,-0.093163580,0.235490034,0.016590287,-0.048400917, 
 0.164004024,0.009760047,0.177393526,-0.051797659,0.081804319,0.169090244)
 
names(active.k.vals) <- c("ENSG00000129993.10", "ENSG00000147883.9", "ENSG00000175806.10", "ENSG00000125249.6", 
 "ENSG00000197747.4", "ENSG00000138744.10", "ENSG00000161800.8", "ENSG00000150403.13", "ENSG00000185803.4", 
 "ENSG00000120262.8", "ENSG00000157600.7", "ENSG00000181751.5")


active.k.vals <- active.k.vals[intersect(names(active.k.vals), rownames(norm.exp))]
risk.score <- RiskEsti(exp.dat=norm.exp, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)

cli.data <- cbind(cli.data, risk.score[rownames(cli.data), ])


# survival plot 
source('/code/Rscript/survival_plot.R')

SurvivalPlot(survival.data=cli.data[, c('Tumor (T) sample ID', 'Overall survial (month)', 'Survial  (1, dead; 0, alive)')], 
 sample.class=cli.data[, c('Tumor (T) sample ID', 'risk.categ')], filename='Gao_Cell_2019_lihc_os.pdf', 
 out.file.path='/result/Section3/')

SurvivalPlot(survival.data=cli.data[, c('Tumor (T) sample ID', 'Recurrence-free survival (month)', 'Recurrence  (1, yes; 0, no)')], 
 sample.class=cli.data[, c('Tumor (T) sample ID', 'risk.categ')], filename='Gao_Cell_2019_lihc_rfs.pdf', 
 out.file.path='/result/Section3/')


saveRDS(cli.data, file='/data/Gao_Cell_2019_risk_score.rds')


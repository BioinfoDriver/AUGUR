
RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
	
	exp.dat <- as.matrix(exp.dat[intersect(rownames(exp.dat), gene.set), ])
	
	risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
	if(!is.null(cut.off)){
		risk.categ <- ifelse(risk.score >= cut.off, 'high risk', 'low risk')
	
	}else{
		risk.categ <- ifelse(risk.score >= median(risk.score), 'high risk', 'low risk')
		
	}
	return(list(risk.score, risk.categ))
}





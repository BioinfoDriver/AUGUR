
# load data
tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')
tcga.lihc.cli.data <- readRDS(file='/data/tcga.lihc.cli.data.rds')
candi.prog.sig.gene <- readRDS(file='/data/candi_prog_sig_gene.rds')


# data integration
rownames(tcga.lihc.vst.exp) <- paste0('ID_', rownames(tcga.lihc.vst.exp))
colnames(tcga.lihc.vst.exp) <- substr(colnames(tcga.lihc.vst.exp), 1, 12)
tcga.lihc.vst.exp <- t(tcga.lihc.vst.exp)

tcga.lihc.cli.data <- tcga.lihc.cli.data[, c('os_time', 'os')]
tcga.lihc.cli.data <- as.matrix(tcga.lihc.cli.data)

candi.prog.sig.gene <- paste0('ID_', candi.prog.sig.gene)


# Regularized Cox Regression
LassoCoxFunction <- function(candi.gene, exp.data, cli.data, n.cv, t.measure=c('deviance', 'C'), file.name){
 library(glmnet)
 library(survival)
 
 t.measure = match.arg(t.measure)
 
 colnames(cli.data) <- c('time', 'status')
 exp.data <- exp.data[, candi.gene]
 
 inters.sams <- intersect(rownames(exp.data), rownames(cli.data))
 exp.data <- exp.data[inters.sams, ]
 cli.data <- cli.data[inters.sams, ]
 
 
 # Does k-fold cross-validation for glmnet,
 set.seed(123)
 cv.fit = cv.glmnet(x = exp.data, y = cli.data, type.measure = t.measure, nfolds=n.cv, family = "cox")
 
 
 # plot the cross-validation curve
 pdf(file.name)
  plot(cv.fit)
 dev.off()

 return(cv.fit)
}

setwd('/result/Section2/')
lasso.cox.res <- LassoCoxFunction(candi.prog.sig.gene, tcga.lihc.vst.exp, 
 tcga.lihc.cli.data, n.cv=10, t.measure='deviance', file.name='cv_curve_deviance.pdf')


# extract non-zero coefficients
est.coef = coef(lasso.cox.res, s = lasso.cox.res$lambda.min)
active.k.vals = est.coef[which(est.coef != 0), ]
names(active.k.vals) <- gsub('ID_', '', names(active.k.vals))

est.coef = coef(lasso.cox.res, s = lasso.cox.res$lambda.1se)
active.k.vals.1se = est.coef[which(est.coef != 0), ]
names(active.k.vals.1se) <- gsub('ID_', '', names(active.k.vals.1se))

save(lasso.cox.res, active.k.vals, active.k.vals.1se, file='/data/lasso_cox_res.RData')


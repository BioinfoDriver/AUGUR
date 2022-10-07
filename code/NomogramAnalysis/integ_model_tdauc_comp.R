
# load data
load(file='/data/four_dataset_integ_model_socre.RData')

# Multiline Time-dependent ROC curve
MultilineTdROC <- function(sur.time, sur.status, risk.scores, t.point, output.file){

 library(timeROC)
 library(survival)
 library(ggplot2)
 library(purrr)
 library(tibble)
 
 ExtractRoc <- function(roc.obj) {
  # extract sens and spec from roc.object as a tibble
  as_tibble(data.frame(FPR=roc.obj$FP[, 2], TPR=roc.obj$TP[, 2]))
 }

 # get roc curves for each predictor
 roc.plot <- risk.scores %>% map(~timeROC(T=sur.time, delta=sur.status, marker = .x, 
 cause=1, times=t.point, ROC=TRUE, iid=TRUE)) %>% map_dfr(ExtractRoc, .id = "pred") %>% 
 ggplot() + geom_abline(slope = 1, intercept = 0, color = "darkgrey") +
 aes(FPR, TPR, color = pred) + geom_line() + theme_bw() + 
 labs(x = "False positive rate", y = "True positive rate")
 
 # save
 ggsave(filename=output.file, plot=roc.plot)
}


setwd('/result/Section6/')

# TCGA
risk.scores <- as.list(tcga.model.score[, c('TNMStage', 'CloBiomarker', 'IntegModel')])
MultilineTdROC(tcga.model.score$os_time, tcga.model.score$os, risk.scores, 365*1, 'tcga_one_year_roc.pdf')
MultilineTdROC(tcga.model.score$os_time, tcga.model.score$os, risk.scores, 365*3, 'tcga_three_year_roc.pdf')
MultilineTdROC(tcga.model.score$os_time, tcga.model.score$os, risk.scores, 365*5, 'tcga_five_year_roc.pdf')

# ICGC
risk.scores <- as.list(icgc.model.score[, c('TNMStage', 'CloBiomarker', 'IntegModel')])
MultilineTdROC(icgc.model.score$os_time, icgc.model.score$os, risk.scores, 365*1, 'icgc_one_year_roc.pdf')
MultilineTdROC(icgc.model.score$os_time, icgc.model.score$os, risk.scores, 365*3, 'icgc_three_year_roc.pdf')

# Gao et al
risk.scores <- as.list(gao.etal.model.score[, c('TNMStage', 'CloBiomarker', 'IntegModel')])
MultilineTdROC(gao.etal.model.score$os_time, gao.etal.model.score$os, risk.scores, 365*1, 'gao_one_year_roc.pdf')
MultilineTdROC(gao.etal.model.score$os_time, gao.etal.model.score$os, risk.scores, 365*3, 'gao_three_year_roc.pdf')

# GSE14520
risk.scores <- as.list(lci.model.score[, c('TNMStage', 'CloBiomarker', 'IntegModel')])
MultilineTdROC(lci.model.score$os_time, lci.model.score$os, risk.scores, 365*1, 'lci_one_year_roc.pdf')
MultilineTdROC(lci.model.score$os_time, lci.model.score$os, risk.scores, 365*3, 'lci_three_year_roc.pdf')
MultilineTdROC(lci.model.score$os_time, lci.model.score$os, risk.scores, 365*5, 'lci_five_year_roc.pdf')



# Time-dependent ROC curve analysis
TimeDependentROC <- function(clinical.data, sur.time, sur.status, models, t.points){

 library(timeROC)
 library(survival)

 time.roc.auc <- lapply(models, function(x){
  
  time.roc <- timeROC(T = sur.time, delta = sur.status, marker = clinical.data[, x], cause=1,
   times = t.points, ROC = TRUE, iid = TRUE)
  
  # Area under the curve
  auc <- time.roc$AUC
  # Confidence intervals
  auc.ci95 <- confint(time.roc, level = 0.95)$CI_AUC/100
 
  roc.auc <- data.frame(cbind(auc, auc.ci95))
  roc.auc$auc_ci <- paste0(sprintf('%0.2f', roc.auc$auc), 
  '(', sprintf('%0.2f', roc.auc$X2.5.), '-', sprintf('%0.2f', roc.auc$X97.5.), ')')
  
  roc.auc$model <- x
  return(roc.auc)
 })
 
 return(time.roc.auc)
}

risk.sig <- c('TNMStage', 'CloBiomarker', 'IntegModel')
# TCGA
TimeDependentROC(tcga.model.score, tcga.model.score$os_time, tcga.model.score$os, 
 risk.sig, t.points = c(1*365, 3*365, 5*365))

# ICGC
TimeDependentROC(icgc.model.score, icgc.model.score$os_time, icgc.model.score$os, 
 risk.sig, t.points = c(1*365, 3*365))

# Gao et al
TimeDependentROC(gao.etal.model.score, gao.etal.model.score$os_time, gao.etal.model.score$os, 
 risk.sig, t.points = c(1*365, 3*365))

# GSE14520
TimeDependentROC(lci.model.score, lci.model.score$os_time, lci.model.score$os, 
 risk.sig, t.points = c(1*365, 3*365, 5*365))



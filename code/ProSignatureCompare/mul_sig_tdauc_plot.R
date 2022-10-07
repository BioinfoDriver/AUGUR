

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

clinical.data <- readRDS(file='/data/tcga_sig_risk_score.rds')

risk.scores <- as.list(clinical.data[, c('risk.score', 'PMID_32198063', 'PMID_30885723', 
 'PMID_35123387', 'PMID_31335995', 'PMID_31695766','PMID_33033585', 'PMID_34676211', 'PMID_34277622',
 'PMID_33777771', 'PMID_34975331', 'PMID_25666192','PMID_22105560')])

setwd('/result/Section5/')
MultilineTdROC(clinical.data$os_time, clinical.data$os, risk.scores, 365*1, 'one_year_roc_comp.pdf')
MultilineTdROC(clinical.data$os_time, clinical.data$os, risk.scores, 365*3, 'three_year_roc_comp.pdf')
MultilineTdROC(clinical.data$os_time, clinical.data$os, risk.scores, 365*5, 'five_year_roc_comp.pdf')




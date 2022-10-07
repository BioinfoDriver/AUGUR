

# Time-dependent ROC curve analysis
TimeDependentROC <- function(sur.time, sur.status, risk.score, t.points, if.plot=TRUE, time.labels, plot.file.name){

 library(timeROC)
 library(survival)
 library(ggplot2)

 time.roc <- timeROC(T = sur.time, delta = sur.status, marker = risk.score, cause=1,
  times = t.points, ROC = TRUE, iid = TRUE)
 
 # Area under the curve
 auc <- time.roc$AUC
 # Confidence intervals
 auc.ic95 <- confint(time.roc, level = 0.95)$CI_AUC


 # Plotting time-dependent ROC curves
 if(if.plot){
  
  time.roc.sta <- data.frame(TP_1 = time.roc$TP[, 1], FP_1 = time.roc$FP[, 1],
   TP_2 = time.roc$TP[, 2], FP_2 = time.roc$FP[, 2], TP_3 = time.roc$TP[, 3], FP_3 = time.roc$FP[, 3])

  # plot
  roc.plot <- ggplot(data = time.roc.sta) +
   geom_line(aes(x = FP_1, y = TP_1), size = 1, color = "#BC3C29FF") +
   geom_line(aes(x = FP_2, y = TP_2), size = 1, color = "#0072B5FF") +
   geom_line(aes(x = FP_3, y = TP_3), size = 1, color = "#E18727FF") +
   geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
   theme_bw() +
   annotate("text", x = 0.75, y = 0.15, size = 3.5, color = "#BC3C29FF",
    label = paste0("AUC at ", time.labels[1], ":", sprintf("%.3f", auc[1]), 
     " (", sprintf("%.3f", auc.ic95[1, 1]), "-", sprintf("%.3f", auc.ic95[1, 2]), ")")
   ) +
   annotate("text", x = 0.75, y = 0.10, size = 3.5, color = "#0072B5FF",
    label = paste0("AUC at ", time.labels[2], ":", sprintf("%.3f", auc[2]), 
     " (", sprintf("%.3f", auc.ic95[2, 1]), "-", sprintf("%.3f", auc.ic95[2, 2]), ")")
   ) +
   annotate("text", x = 0.75, y = 0.05, size = 3.5, color = "#E18727FF",
    label = paste0("AUC at ", time.labels[3], ":", sprintf("%.3f", auc[3]), 
     " (", sprintf("%.3f", auc.ic95[3, 1]), "-", sprintf("%.3f", auc.ic95[3, 2]), ")")
   ) +
   labs(x = "False positive rate", y = "True positive rate") +
   theme(
     axis.text = element_text(face = "bold", size = 11, color = "black"),
     axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
     axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
   )
  
  ggsave(filename=plot.file.name, plot=roc.plot)
  }else{
  
   return(list('auc'=auc, 'auc.ic95'=auc.ic95))
  }
}


# Time-dependent ROC curve analysis
# TCGA
tcga.risk.score <- readRDS(file='/data/tcga_risk_score.rds')
TimeDependentROC(sur.time=tcga.risk.score$os_time, sur.status=tcga.risk.score$os, risk.score=tcga.risk.score$risk.score, 
 t.points=c(1*365, 3*365, 5*365), if.plot=TRUE, time.labels=c('1 years', '3 years','5 years'), 
 plot.file.name='/result/Section3/tcga_os_timeROC.pdf')


TimeDependentROC(sur.time=tcga.risk.score$pfi_time, sur.status=tcga.risk.score$pfi, risk.score=tcga.risk.score$risk.score, 
 t.points=c(1*365, 3*365, 5*365), if.plot=TRUE, time.labels=c('1 years', '3 years','5 years'), 
 plot.file.name='/result/Section3/tcga_pfi_timeROC.pdf')


tcga.risk.dss <- subset(tcga.risk.score, !is.na(dss) & !is.na(dss_time))
TimeDependentROC(sur.time=tcga.risk.dss$dss_time, sur.status=tcga.risk.dss$dss, risk.score=tcga.risk.dss$risk.score, 
 t.points=c(1*365, 3*365, 5*365), if.plot=TRUE, time.labels=c('1 years', '3 years','5 years'), 
 plot.file.name='/result/Section3/tcga_dss_timeROC.pdf')


tcga.risk.dfi <- subset(tcga.risk.score, !is.na(dfi) & !is.na(dfi_time))
TimeDependentROC(sur.time=tcga.risk.dfi$dfi_time, sur.status=tcga.risk.dfi$dfi, risk.score=tcga.risk.dfi$risk.score, 
 t.points=c(1*365, 3*365, 5*365), if.plot=TRUE, time.labels=c('1 years', '3 years','5 years'), 
 plot.file.name='/result/Section3/tcga_dfi_timeROC.pdf')


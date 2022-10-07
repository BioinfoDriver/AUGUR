


RiskScoreDescPlot <- function(risk.score, plot.file.name){
 
 library(ggpubr)
 colnames(risk.score) <- c('patient.id', 'risk.score', 'sur.status')
 risk.score <- subset(risk.score, !is.na(sur.status))
 risk.score$sur.status <- factor(risk.score$sur.status)
 
 risk.plot <- ggbarplot(data=risk.score, x='patient.id', y='risk.score', color=NA, fill='sur.status', sort.by.groups = FALSE, 
  xlab='Risk score for every patient', ylab='Risk score of the CEB-based classifier', sort.val='desc') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  guides(fill=guide_legend(title=NULL))

 ggsave(filename=plot.file.name, plot=risk.plot)
}


tcga.risk.score <- readRDS(file='/data/tcga_risk_score.rds')
icgc.risk.score <- readRDS(file='/data/icgc_risk_score.rds')
GSE144269.risk.score <- readRDS(file='/data/GSE144269_risk_score.rds')
gao.cell.2019.risk.score <- readRDS(file='/data/Gao_Cell_2019_risk_score.rds')
lci.risk.score <- readRDS(file='/data/lci_risk_score.rds')
E_TABM_36.risk.score <- readRDS(file='/data/E_TABM_36_risk_score.rds')
lec.cli.data <- readRDS(file='/data/lec_risk_score.rds')



setwd('/result/Section3/')
RiskScoreDescPlot(tcga.risk.score[, c('patient_id', 'risk.score', 'os')], 'tcga_os_risk.pdf')
RiskScoreDescPlot(tcga.risk.score[, c('patient_id', 'risk.score', 'dss')], 'tcga_dss_risk.pdf')
RiskScoreDescPlot(tcga.risk.score[, c('patient_id', 'risk.score', 'pfi')], 'tcga_pfi_risk.pdf')
RiskScoreDescPlot(tcga.risk.score[, c('patient_id', 'risk.score', 'dfi')], 'tcga_dfi_risk.pdf')

RiskScoreDescPlot(icgc.risk.score[, c('icgc_donor_id', 'risk.score', 'donor_vital_status')], 'icgc_os_risk.pdf')
RiskScoreDescPlot(GSE144269.risk.score[, c('LHC_ID', 'risk.score', 'survival.status')], 'GSE144269_os_risk.pdf')


RiskScoreDescPlot(gao.cell.2019.risk.score[, c('Tumor (T) sample ID', 'risk.score', 'Survial  (1, dead; 0, alive)')], 'gao.cell.2019_os_risk.pdf')
RiskScoreDescPlot(gao.cell.2019.risk.score[, c('Tumor (T) sample ID', 'risk.score', 'Recurrence  (1, yes; 0, no)')], 'gao.cell.2019_rfs_risk.pdf')

RiskScoreDescPlot(lci.risk.score[, c('Affy_GSM', 'risk.score', 'Survival.status')], 'lci_os_risk.pdf')
RiskScoreDescPlot(lci.risk.score[, c('Affy_GSM', 'risk.score', 'Recurr.status')], 'lci_rfs_risk.pdf')

RiskScoreDescPlot(E_TABM_36.risk.score[, c('PatientID', 'risk.score', 'os_status')], 'E_TABM_36_os_risk.pdf')
RiskScoreDescPlot(E_TABM_36.risk.score[, c('PatientID', 'risk.score', 'pfs_status')], 'E_TABM_36_rfs_risk.pdf')

RiskScoreDescPlot(lec.cli.data[, c('Array', 'risk.score', 'OS_Status')], 'lec_os_risk.pdf')


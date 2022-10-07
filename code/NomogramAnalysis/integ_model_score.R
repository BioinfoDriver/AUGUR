

CoxlinearPredictor <- function(cli.data, models){
 
 library(survival)
 risk.scores <- lapply(models, function(model){
 
 ## Fit a Cox model
 coxph <- coxph(formula = as.formula(paste('Surv(os_time, os)~', model)), data = cli.data)
					
 ## Obtain the linear predictor
 risk.score <- predict(coxph, type = "lp")
 return(risk.score)
 })

 # Merge
 risk.scores <- do.call(cbind, risk.scores)
 colnames(risk.scores) <- names(models)
 
 cli.data <- cbind(cli.data, risk.scores)
 return(cli.data)
}


##################TCGA
# load data
tcga.lihc.cli.char <- readRDS(file='/data/tcga.lihc.cli.char.rds')
tcga.lihc.risk.score <- readRDS(file='/data/tcga_risk_score.rds')

# data prepare
sig.score <- tcga.lihc.risk.score[, c('patient_id', 'os', 'os_time', 'risk.score')]
cli.data <- tcga.lihc.cli.char[, c('bcr_patient_barcode', 'pathologic_stage_category')]
cli.data <- subset(cli.data, !is.na(pathologic_stage_category))
cli.data <- merge(sig.score, cli.data, by.x='patient_id', by.y='bcr_patient_barcode')


models <- list(TNMStage='pathologic_stage_category', CloBiomarker='risk.score', 
 IntegModel='pathologic_stage_category + risk.score')
tcga.model.score <- CoxlinearPredictor(cli.data, models)



##################ICGC
# load data
icgc.lihc.risk.score <- readRDS(file='/data/icgc_risk_score.rds')

# data prepare
cli.data <- icgc.lihc.risk.score[, c('icgc_donor_id', 'donor_vital_status', 
 'donor_survival_time', 'risk.score', 'pathologic_stage_category')]
colnames(cli.data) <- c('patient_id', 'os', 'os_time', 'risk.score', 'pathologic_stage')


models <- list(TNMStage='pathologic_stage', CloBiomarker='risk.score', 
 IntegModel='pathologic_stage + risk.score')
icgc.model.score <- CoxlinearPredictor(cli.data, models)


################## Gao et al.
# load data
lihc.risk.score <- readRDS(file='/data/Gao_Cell_2019_risk_score.rds')

# data prepare
cli.data <- lihc.risk.score[, c("Tumor (T) sample ID", "Survial  (1, dead; 0, alive)", 
 "Overall survial (month)", "risk.score", "TNM stage")]
colnames(cli.data) <- c('patient_id', 'os', 'os_time', 'risk.score', 'pathologic_stage')

cli.data$pathologic_stage <- ifelse(cli.data$pathologic_stage %in% c('IA', 'IB', 'II'), 'I/II', 'III/IV')


models <- list(TNMStage='pathologic_stage', CloBiomarker='risk.score', 
 IntegModel='pathologic_stage + risk.score')
gao.etal.model.score <- CoxlinearPredictor(cli.data, models)



#################LCI
# load data
lci.risk.score <- readRDS(file='/data/lci_risk_score.rds')

# data prepare
cli.data <- lci.risk.score[, c("LCS.ID", "Survival.status", "Survival.months", "risk.score", "TNM.staging")]
colnames(cli.data) <- c('patient_id', 'os', 'os_time', 'risk.score', 'pathologic_stage')

cli.data <- subset(cli.data, pathologic_stage != '.')
cli.data$pathologic_stage <- ifelse(cli.data$pathologic_stage %in% c('I', 'II'), 'I/II', 'III/IV')


models <- list(TNMStage='pathologic_stage', CloBiomarker='risk.score', 
 IntegModel='pathologic_stage + risk.score')
lci.model.score <- CoxlinearPredictor(cli.data, models)


# save
save(tcga.model.score, icgc.model.score, gao.etal.model.score, lci.model.score, 
 file='/data/four_dataset_integ_model_socre.RData')


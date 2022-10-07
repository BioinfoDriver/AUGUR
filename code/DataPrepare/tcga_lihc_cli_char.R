

setwd('/data/OriginalData/PanCanAtlas')
cli.data <- read.csv(file='clinical_PANCAN_patient_with_followup.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)

# LIHC clinical characteristics
lihc.cli.data <- subset(cli.data, acronym == 'LIHC')

index.null <- sapply(lihc.cli.data, function(x) all(is.na(x)) | all(x == '') | all(x == '[Not Available]') | all(x=='[Not Applicable]'))
lihc.cli.data <- lihc.cli.data[, !index.null]

# select characteristic
lihc.cli.data <- lihc.cli.data[, c(
"bcr_patient_barcode",
"age_at_initial_pathologic_diagnosis",
"gender",
"race",
"ethnicity",
"height",
"weight",

"histological_type",
"pathologic_T",
"pathologic_M",
"pathologic_N",
"pathologic_stage",
"system_version",

"neoplasm_histologic_grade",
"residual_tumor",
"child_pugh_classification_grade",
"eastern_cancer_oncology_group",
"vascular_tumor_cell_invasion_type", 

"history_hepato_carcinoma_risk_factor",
"history_hepato_carcinoma_risk_factor_other",


"adjacent_hepatic_tissue_inflammation_extent_type",
"liver_fibrosis_ishak_score_category",
   
"history_of_neoadjuvant_treatment",
"relative_family_cancer_history_ind_3",
"prior_dx",

"new_tumor_event_after_initial_treatment",
"person_neoplasm_cancer_status",
"postoperative_rx_tx",

"laboratory_procedure_prothrombin_time_result_value",
"laboratory_procedure_alpha_fetoprotein_outcome_value",
"laboratory_procedure_albumin_result_specified_value",
"lab_procedure_platelet_result_specified_value",
"hematology_serum_creatinine_laboratory_result_value_in_mg_dl"                             
)]

# replace 
lihc.cli.data[lihc.cli.data=='[Not Available]'] <- NA

# converting multiple columns from character to numeric format 
library(dplyr)
lihc.cli.data <- lihc.cli.data %>% mutate_at(c(2, 6, 7, 29:33), as.numeric)


# disposal characteristic 
lihc.cli.data <- lihc.cli.data[, c(
"bcr_patient_barcode",
"age_at_initial_pathologic_diagnosis",
"gender",
"race",
"height",
"weight",
"pathologic_stage",
"neoplasm_histologic_grade",
"residual_tumor",
"child_pugh_classification_grade",
"vascular_tumor_cell_invasion_type", 
"adjacent_hepatic_tissue_inflammation_extent_type",
"liver_fibrosis_ishak_score_category",
"relative_family_cancer_history_ind_3",
"laboratory_procedure_alpha_fetoprotein_outcome_value",
'histological_type',
'history_of_neoadjuvant_treatment',
'prior_dx')]

colnames(lihc.cli.data) <- c("bcr_patient_barcode","age","gender","race","height","weight","pathologic_stage",
"histologic_grade","residual_tumor","child_pugh_grade","vascular_tumor_invasion", "inflammation_extent","ishak_score",
"family_cancer_history","alpha_fetoprotein",'histological_type','history_of_neoadjuvant_treatment','prior_dx')

lihc.cli.data$age_category <- ifelse(lihc.cli.data$age>=60, '≥60', '<60')
lihc.cli.data$age_category <- factor(lihc.cli.data$age_category, levels=c('<60', '≥60'))

lihc.cli.data$gender_category <- factor(lihc.cli.data$gender, levels=c('FEMALE', 'MALE'))

lihc.cli.data$race_category <- lihc.cli.data$race
lihc.cli.data$race_category[lihc.cli.data$race_category %in% c('[Unknown]', '[Not Evaluated]')] <- NA
lihc.cli.data$race_category <- ifelse(lihc.cli.data$race_category == 'ASIAN', 'Asian', 'Not Asian')
lihc.cli.data$race_category <- factor(lihc.cli.data$race_category, levels=c('Not Asian', 'Asian'))

lihc.cli.data$BMI <- lihc.cli.data$weight/((lihc.cli.data$height/100)^2)
lihc.cli.data$BMI_category <- ifelse(lihc.cli.data$BMI>=25, '≥25', '<25')
lihc.cli.data$BMI_category <- factor(lihc.cli.data$BMI_category, levels=c('<25', '≥25'))

lihc.cli.data$pathologic_stage_category <- ifelse(lihc.cli.data$pathologic_stage == 'Stage I' |
 lihc.cli.data$pathologic_stage == 'Stage II', 'I/II', 'III/IV')
lihc.cli.data$pathologic_stage_category <- factor(lihc.cli.data$pathologic_stage_category, levels=c('I/II', 'III/IV'))

lihc.cli.data$histologic_grade_category <- ifelse(lihc.cli.data$histologic_grade == 'G1' |
 lihc.cli.data$histologic_grade == 'G2', '1/2', '3/4')
lihc.cli.data$histologic_grade_category <- factor(lihc.cli.data$histologic_grade_category, levels=c('1/2', '3/4'))

lihc.cli.data$residual_tumor_category <- lihc.cli.data$residual_tumor
lihc.cli.data$residual_tumor_category[lihc.cli.data$residual_tumor_category == 'RX'] <- NA
lihc.cli.data$residual_tumor_category <- ifelse(lihc.cli.data$residual_tumor_category == 'R0', 'R0', 'R1/2')
lihc.cli.data$residual_tumor_category <- factor(lihc.cli.data$residual_tumor_category, levels=c('R0', 'R1/2'))

lihc.cli.data$child_pugh_grade_category <- lihc.cli.data$child_pugh_grade
lihc.cli.data$child_pugh_grade_category[lihc.cli.data$child_pugh_grade_category == '[Unknown]'] <- NA
lihc.cli.data$child_pugh_grade_category <- ifelse(lihc.cli.data$child_pugh_grade_category == 'A', 'A', 'B/C')
lihc.cli.data$child_pugh_grade_category <- factor(lihc.cli.data$child_pugh_grade_category, levels=c('A', 'B/C'))

lihc.cli.data$vascular_tumor_invasion_category <- lihc.cli.data$vascular_tumor_invasion
lihc.cli.data$vascular_tumor_invasion_category[lihc.cli.data$vascular_tumor_invasion_category == '[Unknown]'] <- NA
lihc.cli.data$vascular_tumor_invasion_category <- 
ifelse(lihc.cli.data$vascular_tumor_invasion_category == 'None', 'None', 'Micro/Macro')
lihc.cli.data$vascular_tumor_invasion_category <- 
factor(lihc.cli.data$vascular_tumor_invasion_category, levels=c('None', 'Micro/Macro'))

lihc.cli.data$inflammation_extent_category <- lihc.cli.data$inflammation_extent
lihc.cli.data$inflammation_extent_category[lihc.cli.data$inflammation_extent_category == '[Unknown]'] <- NA
lihc.cli.data$inflammation_extent_category <- ifelse(lihc.cli.data$inflammation_extent_category == 'None', 'None', 'Mild/Severe')
lihc.cli.data$inflammation_extent_category <- 
factor(lihc.cli.data$inflammation_extent_category, levels=c('None', 'Mild/Severe'))

lihc.cli.data$ishak_score_category <- lihc.cli.data$ishak_score
lihc.cli.data$ishak_score_category[lihc.cli.data$ishak_score_category == '[Unknown]'] <- NA
lihc.cli.data$ishak_score_category <- ifelse(lihc.cli.data$ishak_score_category == "6 - Established Cirrhosis" |
 lihc.cli.data$ishak_score_category == "5 - Nodular Formation and Incomplete Cirrhosis", '5-6', '0-4')
lihc.cli.data$ishak_score_category <- factor(lihc.cli.data$ishak_score_category, levels=c('0-4', '5-6'))

lihc.cli.data$family_cancer_history_category <- lihc.cli.data$family_cancer_history
lihc.cli.data$family_cancer_history_category[lihc.cli.data$family_cancer_history_category == '[Unknown]'] <- NA
lihc.cli.data$family_cancer_history_category <- factor(lihc.cli.data$family_cancer_history_category, levels=c('No', 'Yes'))

lihc.cli.data$alpha_fetoprotein <- ifelse(lihc.cli.data$alpha_fetoprotein>=300, '≥300', '<300')
lihc.cli.data$alpha_fetoprotein <- factor(lihc.cli.data$alpha_fetoprotein, levels=c('<300', '≥300'))

saveRDS(lihc.cli.data, file='/data/tcga.lihc.cli.char.rds')

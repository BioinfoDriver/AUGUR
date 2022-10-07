
# clinical data
setwd('/data/OriginalData/ICGC_LIRI_JP')
donor <- read.csv(file='donor.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)

library(dplyr)
donor <- donor %>% select(!c('project_code', 'study_donor_involved_in','disease_status_last_followup',
 'donor_relapse_type','donor_age_at_enrollment', 'donor_age_at_last_followup', 'donor_relapse_interval',
 'donor_diagnosis_icd10','donor_tumour_stage_at_diagnosis_supplemental', 'donor_interval_of_last_followup'))

# merge
icgc.rnaseq.sample <- readRDS(file='/data/icgc_rnaseq_filter_sample.rds')
icgc.rnaseq.sample <- subset(icgc.rnaseq.sample, specimen_type == 'Tumor')

icgc.lihc.cli.data <- merge(donor, icgc.rnaseq.sample, by = 'icgc_donor_id', all = TRUE)
icgc.lihc.cli.data$filter <- ifelse(is.na(icgc.lihc.cli.data$specimen_type), TRUE, FALSE)


# Clinical and pathological information of the 300 liver cancers
cli.patho <- read.csv(file='ClinicalAndPathologicalInformation.txt', header=TRUE, sep='\t')

# curated
cli.patho <- subset(cli.patho, ID %in% c(icgc.lihc.cli.data$submitted_donor_id.x, 
 'RK006_1', 'RK180_1', 'RK233_1', 'RK261_1', 'RK297_1', 'RK338_1'))
cli.patho$ID <- substr(cli.patho$ID, 1, 5)

icgc.lihc.cli.data <- merge(icgc.lihc.cli.data, cli.patho, by.x = 'submitted_donor_id.x', by.y = 'ID', all.x = TRUE)

# 两个数据来源，不一致的样本有, 因此以ICGC数据为参考
# Sex: RK105, RK106, RK135
# Age: RK105, RK106
# Stage: RK019, RK105, RK106
# Grade: RK004, RK180, "RK024", "RK152", "RK162", "RK163", "RK175", "RK285", RK148, RK243, RK280

# disposal characteristic
icgc.lihc.cli.data$age_category <- ifelse(icgc.lihc.cli.data$donor_age_at_diagnosis>=60, '≥60', '<60')
icgc.lihc.cli.data$age_category <- factor(icgc.lihc.cli.data$age_category, levels=c('<60', '≥60'))

icgc.lihc.cli.data$gender_category <- factor(icgc.lihc.cli.data$donor_sex, levels=c('female', 'male'))

icgc.lihc.cli.data$pathologic_stage_category <- 
 ifelse(icgc.lihc.cli.data$donor_tumour_stage_at_diagnosis %in% c(1, 2), 'I/II', 'III/IV')
icgc.lihc.cli.data$pathologic_stage_category <- factor(icgc.lihc.cli.data$pathologic_stage_category, levels=c('I/II', 'III/IV'))


icgc.lihc.cli.data$histologic_grade_category <- 
 ifelse(icgc.lihc.cli.data$tumour_grade %in% c('I', 'II', 'I-II', 'II-I'), '1/2', '3/4')
icgc.lihc.cli.data$histologic_grade_category <- factor(icgc.lihc.cli.data$histologic_grade_category, levels=c('1/2', '3/4'))

icgc.lihc.cli.data$tumor_size_category <- ifelse(icgc.lihc.cli.data$Tumor.size..mm. >= 30, '≥30', '<30')
icgc.lihc.cli.data$tumor_size_category <- factor(icgc.lihc.cli.data$tumor_size_category, levels=c('<30', '≥30'))

icgc.lihc.cli.data$portal_invasion_category <- ifelse(icgc.lihc.cli.data$Portal.vein.invasion == 0, 'No', 'Yes')
icgc.lihc.cli.data$portal_invasion_category <- factor(icgc.lihc.cli.data$portal_invasion_category, levels=c('No', 'Yes'))

icgc.lihc.cli.data$portal_invasion_category <- ifelse(icgc.lihc.cli.data$Portal.vein.invasion == 0, 'No', 'Yes')
icgc.lihc.cli.data$portal_invasion_category <- factor(icgc.lihc.cli.data$portal_invasion_category, levels=c('No', 'Yes'))

icgc.lihc.cli.data$hepatic_invasion_category <- ifelse(icgc.lihc.cli.data$Hepatic.vein.invasion == 0, 'No', 'Yes')
icgc.lihc.cli.data$hepatic_invasion_category <- factor(icgc.lihc.cli.data$hepatic_invasion_category, levels=c('No', 'Yes'))

icgc.lihc.cli.data$bile_invasion_category <- ifelse(icgc.lihc.cli.data$Bile.duct.invasion == 0, 'No', 'Yes')
icgc.lihc.cli.data$bile_invasion_category <- factor(icgc.lihc.cli.data$bile_invasion_category, levels=c('No', 'Yes'))

icgc.lihc.cli.data$fibrosisc_category <- ifelse(icgc.lihc.cli.data$Liver.fibrosisc == 4, 'Yes', 'No')
icgc.lihc.cli.data$fibrosisc_category <- factor(icgc.lihc.cli.data$fibrosisc_category, levels=c('No', 'Yes'))

icgc.lihc.cli.data$Virus.infection <- trimws(icgc.lihc.cli.data$Virus.infection)
icgc.lihc.cli.data$virus_infection_category <- ifelse(icgc.lihc.cli.data$Virus.infection=='NBNC', 'NBNC', 'HBV/HCV')
icgc.lihc.cli.data$virus_infection_category <- factor(icgc.lihc.cli.data$virus_infection_category, levels=c('NBNC', 'HBV/HCV'))

icgc.lihc.cli.data$alcohol_intake_category <- ifelse(icgc.lihc.cli.data$Alcohol.intake==0, 'No', 'Yes')
icgc.lihc.cli.data$alcohol_intake_category <- factor(icgc.lihc.cli.data$alcohol_intake_category, levels=c('No', 'Yes'))

saveRDS(icgc.lihc.cli.data, file='/data/icgc_linc_cli_data.rds')





library(tibble)
library(data.table)
library(tidyr)
library(dplyr)
library(bit64)

# RNA-seq
setwd('/data/ICGC_LIRI_JP')
icgc.lihc.exp <- fread(file='exp_seq.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
icgc.lihc.exp <- mutate(icgc.lihc.exp, across(where(is.integer64), as.numeric))

icgc.lihc.fpkm <- spread(icgc.lihc.exp[, c('icgc_sample_id', 'gene_id', 'normalized_read_count')], 
 gene_id, normalized_read_count)
icgc.lihc.count <- spread(icgc.lihc.exp[, c('icgc_sample_id', 'gene_id', 'raw_read_count')], 
 gene_id, raw_read_count)

icgc.lihc.fpkm <- icgc.lihc.fpkm %>% column_to_rownames(var = "icgc_sample_id") 
icgc.lihc.count <- icgc.lihc.count %>% column_to_rownames(var = "icgc_sample_id") 


icgc.lihc.fpkm <- as.data.frame(t(icgc.lihc.fpkm))
icgc.lihc.count <- as.data.frame(t(icgc.lihc.count))

# sample annotation of RNA-seq data (icgc_donor_id:232, icgc_sample_id:445)
icgc.rnaseq.sample <- icgc.lihc.exp[, c('icgc_specimen_id', 'icgc_sample_id', 'submitted_sample_id', 
 'analysis_id', 'total_read_count')]
icgc.rnaseq.sample <- icgc.rnaseq.sample[!duplicated(icgc.rnaseq.sample), ]


# clinical data
specimen <- read.csv(file='specimen.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)
icgc.rnaseq.sample <- merge(icgc.rnaseq.sample, specimen, by = 'icgc_specimen_id')
# 243 liver cancer samples from 232 patients
# > table(icgc.rnaseq.sample[, c('specimen_type', 'tumour_histological_type')])
                                           # tumour_histological_type
# specimen_type                                   8130/3 8160/3 8170/3 8180/3
  # Metastatic tumour - additional metastatic   0      3      0      0      0
  # Normal - solid tissue                     197      0      0      0      0
  # Normal - tissue adjacent to primary         5      0      0      0      0
  # Primary tumour - solid tissue               0      3     20    210      7

# including three with metastatic tumor, 30 with ICC or cHCC/CC
# ICD-O International Classification of Diseases for Oncology
# https://www.genome.jp/kegg-bin/get_htext#A5
# 8130/3  Papillary transitional cell carcinoma
# 8160/3  Cholangiocarcinoma
# 8170/3  Hepatocellular carcinoma, NOS
# 8180/3  Combined hepatocellular carcinoma and cholangiocarcinoma
icgc.rnaseq.sample <- subset(icgc.rnaseq.sample, specimen_type != 'Metastatic tumour - additional metastatic')
icgc.rnaseq.sample <- subset(icgc.rnaseq.sample, !tumour_histological_type %in% c('8130/3', '8160/3', '8180/3'))


# 7 duplicated samples
# remove 4 deplicated samples(样本编号不同, 测序总read数目和基因表达值一模一样)
# RK261_Cancer2(SA595567) == RK261_Cancer(SA595566) (DO50814)
# RK233_Cancer2(SA595535) == RK233_Cancer(SA595534) (DO50793)
# RK006_Cancer2(SA594231) == RK006_Cancer(SA594223) (DO45096)
# RK180_Cancer2(SA595463) == RK180_Cancer(SA595462) (DO45273)

# remove 3 duplicated samples with low tumor cell percentage
# RK189_Cancer(SA595474)/RK036_Cancer(SA560741) (DO45141)
# RK240_Cancer(SA595543)/RK179_Cancer(SA595460)/RK035_Cancer(SA594432) (DO45139)


icgc.rnaseq.sample <- subset(icgc.rnaseq.sample, !(icgc_sample_id %in% 
 c('SA595567', 'SA595535', 'SA594231', 'SA595463', 'SA595474', 'SA595543', 'SA595460')))
# 删除后患者数目为(icgc_donor_id:226, icgc_sample_id:405, tumor:203, normal:202, 其中176个患者具有配对正常)

icgc.lihc.fpkm <- dplyr::select(icgc.lihc.fpkm, icgc.rnaseq.sample$icgc_sample_id)
icgc.lihc.count <- dplyr::select(icgc.lihc.count, icgc.rnaseq.sample$icgc_sample_id)


icgc.rnaseq.sample <- icgc.rnaseq.sample[, -c('project_code', 'study_specimen_involved_in', 'specimen_type_other', 
 'specimen_interval', 'specimen_donor_treatment_type', 'specimen_donor_treatment_type_other', 
 'specimen_processing', 'specimen_processing_other', 'specimen_storage', 'specimen_storage_other', 
 'tumour_confirmed', 'specimen_biobank', 'specimen_biobank_id', 'specimen_available', 
 'tumour_grade_supplemental', 'tumour_stage_supplemental', 'digital_image_of_stained_section', 'level_of_cellularity')]

icgc.rnaseq.sample$specimen_type <- 
 ifelse(icgc.rnaseq.sample$specimen_type == 'Primary tumour - solid tissue', 'Tumor', 'Normal')


icgc.lihc.filter.fpkm <- icgc.lihc.fpkm[, icgc.rnaseq.sample$icgc_sample_id]
icgc.lihc.filter.count <- icgc.lihc.count[, icgc.rnaseq.sample$icgc_sample_id]


# save data
saveRDS(icgc.rnaseq.sample, file='/data/icgc_rnaseq_filter_sample.rds')
save(icgc.lihc.fpkm, icgc.lihc.count, file='/data/icgc_lihc_rnaseq.RData')
save(icgc.lihc.filter.fpkm, icgc.lihc.filter.count, file='/data/icgc_lihc_filter_rnaseq.RData')



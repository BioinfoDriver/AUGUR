
# TCGA sample quality annotations file
setwd('/data/OriginalData/PanCanAtlas')
sam.qua.anno <- read.csv(file='merged_sample_quality_annotations.tsv', header=TRUE, sep='\t',stringsAsFactors=FALSE)
sam.qua.anno <- subset(sam.qua.anno, cancer.type == 'LIHC' & platform == 'IlluminaHiSeq_RNASeqV2')# 424

# remove three with recurrent tumor
sam.qua.anno <- subset(sam.qua.anno, substr(aliquot_barcode, 14, 15)!='02')


# Cholangiocarcinoma(2), cholongio-liver mix(2), not representative of a hepatocellular case(1)

# prior Wilms tumor(1), prior prostate adenocarcinoma(3), prior tumor in the back(1)
# prior colon adenocarcinoma(1), prior papillary transitional cell bladder(1), prior bladder and prostate cancer(1)

# Neoadjuvant therapy(2)

# recurrence(1)
sam.qua.anno <- subset(sam.qua.anno, patient_annotation == '')

# Unlikely HCC(1)
sam.qua.anno <- subset(sam.qua.anno, AWG_pathology_exclusion_reason != 'Unlikely HCC')#01(354), 11(44)


# load data
tcga.lihc.cli.data <- readRDS(file='/data/tcga.lihc.cli.data.rds')
tcga.lihc.cli.char <- readRDS(file='/data/tcga.lihc.cli.char.rds')
load(file='/data/tcga_lihc_rnaseq.RData')


tcga.lihc.cli.data <- tcga.lihc.cli.data[intersect(rownames(tcga.lihc.cli.data), sam.qua.anno$patient_barcode), ] 
tcga.lihc.cli.char <- subset(tcga.lihc.cli.char, bcr_patient_barcode %in% sam.qua.anno$patient_barcode)
# history_of_neoadjuvant_treatment, No(354)

# Fibrolamellar Carcinoma[3], Hepatocellular Carcinoma[347], Hepatocholangiocarcinoma (Mixed)[4]
tcga.lihc.cli.data <- subset(tcga.lihc.cli.data, histological_type == 'Hepatocellular Carcinoma')

# # test and validation
tcga.sam.id <- colnames(tcga.lihc.tpm) # 423, 01(371), 02(2), 11(50)
tcga.prim.sam.id <- substr(tcga.sam.id[substr(tcga.sam.id, 14, 15) == '01'], 1, 12)


tcga.lihc.cli.data <- tcga.lihc.cli.data[intersect(rownames(tcga.lihc.cli.data), tcga.prim.sam.id), ] # 347
tcga.lihc.cli.data <- subset(tcga.lihc.cli.data, !is.na(os) & !is.na(os_time)) # 346


# Only patients with a follow-up period longer than 1 month were included for survival analysis. 
# Patients who were analyzed were those who survived longer than 1 month.
# Including only patients with a follow-up of more than 1 month.
# Only patients with a follow-up period longer than 30 days were included. 
tcga.lihc.cli.data <- subset(tcga.lihc.cli.data, os_time >= 30) # 323, delete 23 patients

set.seed(123)
tcga.train.sam.set <- sample(x=rownames(tcga.lihc.cli.data), size=215, replace = FALSE)
tcga.test.sam.set <- setdiff(rownames(tcga.lihc.cli.data), tcga.train.sam.set)

save(tcga.train.sam.set, tcga.test.sam.set, file='/data/tcga_lihc_tv_samples.RData')


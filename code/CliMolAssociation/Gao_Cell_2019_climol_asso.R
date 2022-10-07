
# load data
lihc.risk.score <- readRDS(file='/data/Gao_Cell_2019_risk_score.rds')


# data prepare

# cli.sig.char <- lihc.risk.score[, c('Tumor (T) sample ID', 
 # 'Age', 'Gender', 'Liver cirrhosis (1, yes; 0, no)', 'Tumor number', 'Tumor size (cm)', 'Lymph node metastasis (1, yes; 0, no)', 
 # 'Tumor thrombus (1, yes; 0, no)', 'Tumour enapsulation (1, complete; 0, no)', 'BCLC stage', 'TNM stage', 
 # 'TB, total bilirubin (µmol/L)', 'ALB, albumin, (g/L)', 'ALT, aminoleucine transferase (U/L)', 
 # 'γ-GT, γ-glutamyltransferase (U/L)', 'Preoperative  AFP（ng/mL）', 'risk.categ')]

cli.sig.char <- lihc.risk.score[, c(1, 4, 3, 10:15, 24, 25, 19:23, 45)]
colnames(cli.sig.char) <- c('patientID', 'age', 'gender', 'cirrhosis', 'tumorNumber', 'tumorSize', 'lymphNodeMetastasis',
 'tumorThrombus', 'tumourEnapsulation', 'BCLCstage', 'TNMstage', 'totalBilirubin', 'ALB', 'ALT', 'GGT', 'AFP', 'risk.categ')

cli.sig.char$age <- ifelse(cli.sig.char$age >= 60, '≥60', '<60')
cli.sig.char$tumorNumber <- ifelse(cli.sig.char$tumorNumber >= 2, '≥2', '1')
cli.sig.char$tumorNumber <- factor(cli.sig.char$tumorNumber, levels=c('1', '≥2'))

cli.sig.char$tumorSize <- ifelse(cli.sig.char$tumorSize > 5, '>5', '≤5')
cli.sig.char$tumorSize <- factor(cli.sig.char$tumorSize, levels=c('≤5', '>5'))

cli.sig.char$tumourEnapsulation <- factor(cli.sig.char$tumourEnapsulation, levels=c(1, 0))


cli.sig.char$BCLCstage <- ifelse(cli.sig.char$BCLCstage %in% c('C'), 'C', 'A/B')
cli.sig.char$TNMstage <- ifelse(cli.sig.char$TNMstage %in% c('IA', 'IB', 'II'), 'I/II', 'III/IV')
cli.sig.char$totalBilirubin <- ifelse(cli.sig.char$totalBilirubin >= 20.5, 'Yes', 'No')     
cli.sig.char$ALB <- ifelse(cli.sig.char$ALB < 35, 'Yes', 'No')     
cli.sig.char$ALT <- ifelse(cli.sig.char$ALT > 50, 'Yes', 'No')     
cli.sig.char$GGT <- ifelse(cli.sig.char$GGT >= 40, 'Yes', 'No')     
cli.sig.char$AFP <- ifelse(cli.sig.char$AFP > 300, 'Yes', 'No')     

cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))


# Fisher's Exact Test
p.values <- lapply(colnames(cli.sig.char)[2:16], function(char){
 print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
 print(table(cli.sig.char[, c(char, 'risk.categ')]))
 
 or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
 p.value <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
 return(c(or, p.value))
})

p.values <- as.data.frame(do.call(rbind, p.values))
colnames(p.values) <- c('OR', 'Pvalue')
p.values$Qvalue <- p.adjust(p.values$Pvalue, 'fdr')
rownames(p.values) <- colnames(cli.sig.char)[2:16]


####Plot
library('ggpubr')

afp.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('AFP', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "AFP", color="AFP", ylab='Number of patients', xlab=FALSE, label=TRUE, 
 lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='AFP', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))


ggsave(ggarrange(afp.p, ncol=4, nrow=3, legend = "none"), file='/result/Section4/Gao_Cell_2019_cli_risk_com.pdf')


##########################################Mutation data

setwd('/data/OriginalData/Gao_Cell_2019')
mut.data <- readxl::read_xlsx(path='mmc1.xlsx', sheet = 3, col_names = TRUE)
colnames(mut.data) <- c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'End_Position',
 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Type', 'Function_in_refGene', 'Variant_Classification', 
 'Amino_acid_Change', 'Normal_depth', 'Tumor_depth', 'VAF')


mut.data$Variant_Classification[mut.data$Variant_Classification == '.'] <- 'Splice_Site'

mut.data$Variant_Type[mut.data$Variant_Classification == 'frameshift deletion'] <- 'DEL'
mut.data$Variant_Classification[mut.data$Variant_Classification == 'frameshift deletion'] <- 'Frame_Shift_Del'

mut.data$Variant_Type[mut.data$Variant_Classification == 'frameshift insertion'] <- 'INS'
mut.data$Variant_Classification[mut.data$Variant_Classification == 'frameshift insertion'] <- 'Frame_Shift_Ins'

mut.data$Variant_Type[mut.data$Variant_Classification == 'nonframeshift deletion'] <- 'DEL'
mut.data$Variant_Classification[mut.data$Variant_Classification == 'nonframeshift deletion'] <- 'In_Frame_Del'

mut.data$Variant_Type[mut.data$Variant_Classification == 'nonframeshift insertion'] <- 'INS'
mut.data$Variant_Classification[mut.data$Variant_Classification == 'nonframeshift insertion'] <- 'In_Frame_Ins'

mut.data$Variant_Classification[mut.data$Variant_Classification %in% 
 c('stopgain', 'stoploss')] <- 'Nonsense_Mutation'

mut.data$Variant_Classification[mut.data$Variant_Classification == 'nonsynonymous SNV'] <- 'Missense_Mutation'


mut.data$Variant_Type[mut.data$Variant_Classification %in% 
 c('frameshift substitution', 'nonframeshift substitution')] <- 'SNV'
mut.data$Variant_Classification[mut.data$Variant_Classification %in% 
 c('frameshift substitution', 'nonframeshift substitution')] <- 'Missense_Mutation'
  

mut.data$Variant_Type[mut.data$Variant_Type == 'INDEL' & 
 mut.data$Reference_Allele != '-' & mut.data$Tumor_Seq_Allele2 != '-' ] <- 'SNV'

mut.data$Variant_Type[mut.data$Variant_Type == 'INDEL' & 
 mut.data$Reference_Allele == '-' & mut.data$Tumor_Seq_Allele2 != '-' ] <- 'INS'

mut.data$Variant_Type[mut.data$Variant_Type == 'INDEL' & 
 mut.data$Reference_Allele != '-' & mut.data$Tumor_Seq_Allele2 == '-' ] <- 'DEL'


library('maftools')

colnames(cli.sig.char)[1] <- 'Tumor_Sample_Barcode'
mut.maf = read.maf(maf = mut.data, clinicalData = cli.sig.char, verbose = TRUE)



# Changing colors for variant classifications
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c('Frame_Shift_Del', 'Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit',
 'Frame_Shift_Ins', 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

show_genes <- c('TP53', 'NFE2L2','ACVR2A','KEAP1','ARID2','BAP1','RB1','AXIN1','ARID1A','APOB','ALB','CTNNB1',
'BRD7', 'CDKN1A', 'CDKN2A', 'CREB3L3', 'DHX9', 'EEF1A1', 'GNAS','HIST1H1E', 'IDH1', 'IL6ST', 'KRAS', 'LZTR1', 'NRAS',
'NUP113','PIK3CA', 'RPS6KA3', 'SMARCA4', 'TSC2', 'WHSC1', 'XPO1')

cols = RColorBrewer::brewer.pal(n = 10,name = 'Spectral')
risk.cols <- cols[1:2]
names(risk.cols) <- c('high risk', 'low risk')
stage.cols <- cols[3:4]
names(stage.cols) <- c('I/II', 'III/IV')
age.cols <- cols[5:6]
names(age.cols) <- c('<60', '≥60')
ggt.cols <- cols[7:8]
names(ggt.cols) <- c('Normal', 'Abnormal')
afp.cols <- cols[9:10]
names(afp.cols) <- c('Normal', 'Abnormal')

 
annocolors = list(risk.categ = risk.cols, TNMstage = stage.cols, age = age.cols, GGT = ggt.cols, AFP = afp.cols)

pdf('/result/Section4/Gao_Cell_2019_mut_oncoplot.pdf')
oncoplot(maf = mut.maf, genes = show_genes, colors = vc_cols, # minMut=9,
  clinicalFeatures = c('risk.categ','TNMstage', 'age', 'GGT', 'AFP'),
  annotationColor = annocolors, sortByAnnotation = TRUE)
dev.off()


# Clinical enrichment analysis
fab.ce = clinicalEnrichment(maf = mut.maf, clinicalFeature = 'risk.categ', minMut = 5)
# fab.ce$pairwise_comparision[fdr < 0.05]

pdf('/result/Section4/Gao_Cell_2019_mut_enrich.pdf')
 plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.1, geneFontSize = 0.5, annoFontSize = 0.6)
dev.off()



# Fisher's Exact Test
lihc.mut.maf <- as.data.frame(mut.maf@data)
 
lihc.mut.maf <- subset(lihc.mut.maf, !Variant_Classification %in% 
 c("3'Flank", "3'UTR", "5'Flank","5'UTR", "IGR", "Intron", "RNA", "Silent"))

mut.genes <- unique(lihc.mut.maf$Hugo_Symbol)
gene.mut.mat <- sapply(unique(lihc.mut.maf$Tumor_Sample_Barcode), function(patient.id){
  gene.sym <- subset(lihc.mut.maf, Tumor_Sample_Barcode == patient.id)$Hugo_Symbol
  index <- as.numeric(mut.genes %in% gene.sym)
  return(index)
})

rownames(gene.mut.mat) <- mut.genes
colnames(gene.mut.mat) <- unique(lihc.mut.maf$Tumor_Sample_Barcode)
high.fre.mut <- gene.mut.mat[rowSums(gene.mut.mat)/ncol(gene.mut.mat) >= 0.03, ]
high.fre.mut <- as.data.frame(t(high.fre.mut))

cli.sig.char <- subset(cli.sig.char, !is.na(Tumor_Sample_Barcode))
high.fre.mut <- high.fre.mut[cli.sig.char$Tumor_Sample_Barcode, ]

# driver genes with 3% mutation frequence
high.fre.mut <- high.fre.mut[, intersect(show_genes, colnames(high.fre.mut))]
high.fre.mut <- cbind(high.fre.mut, cli.sig.char[, 'risk.categ', FALSE])


p.values <- sapply(setdiff(colnames(high.fre.mut), 'risk.categ'), function(driver.gene){
   
	num <- table(high.fre.mut[, driver.gene])
	char <- names(num)[num >= 5]
	dat <- subset(high.fre.mut, get(driver.gene) %in% char)
	
	if(dplyr::n_distinct(dat[, driver.gene])<=1){
	 p <- NA
	
	}else{
	 stat <- table(dat[, c(driver.gene, 'risk.categ')])
	 p <- fisher.test(stat, alternative = "two.sided")$p.value
	
	}
	
	return(p)

})

      # TP53      KEAP1      ARID2        RB1      AXIN1     ARID1A       APOB        ALB     CTNNB1       BRD7     CDKN2A 
# 0.01504362 0.21025553 0.71770668 0.09846842 0.83574242 0.03634980 0.60674910 0.27811929 0.01480998 1.00000000 0.11443616 
      # GNAS    RPS6KA3    SMARCA4       TSC2 
# 0.02742573 0.71770668 0.27336967 0.05654825

p.adjust(p.values, method='fdr')
     # TP53     KEAP1     ARID2       RB1     AXIN1    ARID1A      APOB       ALB    CTNNB1      BRD7    CDKN2A      GNAS 
# 0.1128271 0.3942291 0.8281231 0.2452203 0.8954383 0.1363117 0.8273851 0.4171789 0.1128271 1.0000000 0.2452203 0.1363117 
  # RPS6KA3   SMARCA4      TSC2 
# 0.8281231 0.4171789 0.1696447

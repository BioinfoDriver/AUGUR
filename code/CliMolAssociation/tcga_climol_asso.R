
##########################################Clinical data
# load data
tcga.lihc.cli.char <- readRDS(file='/data/tcga.lihc.cli.char.rds')
tcga.lihc.risk.score <- readRDS(file='/data/tcga_risk_score.rds')

# data prepare
sig.score <- tcga.lihc.risk.score[, c('patient_id', 'os', 'os_time', 'risk.score', 'risk.categ')]
sig.score$risk.categ <- factor(sig.score$risk.categ, levels=c('low risk', 'high risk'))

cli.char <- tcga.lihc.cli.char[, c('bcr_patient_barcode', 'age_category', 'gender_category', 'race_category', 
'family_cancer_history_category', 'BMI_category', 'pathologic_stage_category', 'histologic_grade_category', 
'residual_tumor_category', 'ishak_score_category', 'child_pugh_grade_category', 'alpha_fetoprotein',
'vascular_tumor_invasion_category','inflammation_extent_category')]


cli.sig.char <- merge(sig.score, cli.char, by.x='patient_id', by.y='bcr_patient_barcode')

# Fisher's Exact Test
p.values <- lapply(colnames(cli.sig.char)[6:18], function(char){
 print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
 print(table(cli.sig.char[, c(char, 'risk.categ')]))
 
 or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
 p.value <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
 return(c(or, p.value))
})

p.values <- as.data.frame(do.call(rbind, p.values))
colnames(p.values) <- c('OR', 'Pvalue')
p.values$Qvalue <- p.adjust(p.values$Pvalue, 'fdr')
rownames(p.values) <- colnames(cli.sig.char)[6:18]


####Plot
library('ggpubr')

stage.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('pathologic_stage_category', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "pathologic_stage_category", color="pathologic_stage_category", 
 ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='Stage', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

grade.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('histologic_grade_category', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "histologic_grade_category", color="histologic_grade_category", 
 ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='Grade', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

vascular.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('vascular_tumor_invasion_category', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "vascular_tumor_invasion_category", color="vascular_tumor_invasion_category", 
 ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='Vascular', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

afp.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('alpha_fetoprotein', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "alpha_fetoprotein", color="alpha_fetoprotein", 
 ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='AFP', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))


ggsave(ggarrange(stage.p, grade.p, vascular.p, afp.p, ncol=4, nrow=3, legend = "none"), 
 file='/result/Section4/tcga_cli_risk_com.pdf')


##########################################Mutation data
pan.mut.maf <- data.table::fread(file='/data/OriginalData/PanCanAtlas/mc3.v0.2.8.PUBLIC.maf', sep='\t')

lihc.mut.maf <- subset(pan.mut.maf, substr(Tumor_Sample_Barcode, 1, 12) %in% cli.sig.char$patient_id)
lihc.mut.maf <- subset(lihc.mut.maf, substr(Tumor_Sample_Barcode, 14, 15) == '01')
lihc.mut.maf <- subset(lihc.mut.maf, FILTER == 'PASS')

lihc.mut.maf$Tumor_Sample_Barcode <- substr(lihc.mut.maf$Tumor_Sample_Barcode, 1, 12)
cli.sig.char <- subset(cli.sig.char, patient_id %in% lihc.mut.maf$Tumor_Sample_Barcode)
colnames(cli.sig.char)[1] <- 'Tumor_Sample_Barcode'

library('maftools')
mut.maf = read.maf(maf = lihc.mut.maf, clinicalData = cli.sig.char, verbose = TRUE)


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
grade.cols <- cols[5:6]
names(grade.cols) <- c('1/2', '3/4')
inva.cols <- cols[7:8]
names(inva.cols) <- c('None', 'Micro/Macro')
afp.cols <- cols[9:10]
names(afp.cols) <- c('<300', '≥300')



annocolors = list(risk.categ = risk.cols, alpha_fetoprotein = afp.cols, pathologic_stage_category = stage.cols, 
 histologic_grade_category = grade.cols, vascular_tumor_invasion_category = inva.cols)

pdf('/result/Section4/tcga_mut_oncoplot.pdf')
oncoplot(maf = mut.maf, genes = show_genes, colors = vc_cols, # minMut=9,
  clinicalFeatures = c('risk.categ','alpha_fetoprotein','pathologic_stage_category', 
  'histologic_grade_category','vascular_tumor_invasion_category'),
  annotationColor = annocolors, sortByAnnotation = TRUE)
dev.off()


# Clinical enrichment analysis
fab.ce = clinicalEnrichment(maf = mut.maf, clinicalFeature = 'risk.categ', minMut = 9)
# fab.ce$pairwise_comparision[fdr < 0.05]

pdf('/result/Section4/tcga_mut_enrich.pdf')
 plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)
dev.off()




# Fisher's Exact Test
lihc.mut.maf <- subset(lihc.mut.maf, !Variant_Classification %in% 
 c("3'Flank", "3'UTR", "5'Flank","5'UTR", "Intron", "RNA", "Silent"))

mut.genes <- unique(lihc.mut.maf$Hugo_Symbol)
gene.mut.mat <- sapply(unique(lihc.mut.maf$Tumor_Sample_Barcode), function(patient.id){
  gene.sym <- subset(lihc.mut.maf, Tumor_Sample_Barcode == patient.id)$Hugo_Symbol
  index <- as.numeric(mut.genes %in% gene.sym)
  return(index)
})

rownames(gene.mut.mat) <- mut.genes
high.fre.mut <- gene.mut.mat[rowSums(gene.mut.mat)/ncol(gene.mut.mat) >= 0.03, ]
high.fre.mut <- as.data.frame(t(high.fre.mut))
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

        # TP53       ACVR2A        KEAP1        ARID2         BAP1          RB1 
# 1.144914e-09 2.183006e-01 5.969958e-01 1.325906e-01 3.187964e-01 5.838000e-03 
       # AXIN1       ARID1A         APOB          ALB       CTNNB1         BRD7 
# 3.852069e-01 1.492890e-01 1.895809e-01 2.937802e-01 2.049415e-03 1.000000e+00 
        # DHX9        LZTR1       PIK3CA      RPS6KA3      SMARCA4         TSC2 
# 1.000000e+00 7.692459e-01 6.085698e-02 7.977789e-01 3.767037e-01 1.015594e-02 

p.adjust(p.values, method='fdr')
        # TP53       ACVR2A        KEAP1        ARID2         BAP1          RB1 
# 2.060845e-08 4.366012e-01 7.675660e-01 3.838861e-01 5.216669e-01 3.502800e-02 
       # AXIN1       ARID1A         APOB          ALB       CTNNB1         BRD7 
# 5.333634e-01 3.838861e-01 4.265571e-01 5.216669e-01 1.844473e-02 1.000000e+00 
        # DHX9        LZTR1       PIK3CA      RPS6KA3      SMARCA4         TSC2 
# 1.000000e+00 8.975013e-01 2.190851e-01 8.975013e-01 5.333634e-01 4.570173e-02 

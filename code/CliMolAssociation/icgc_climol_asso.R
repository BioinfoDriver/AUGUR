##########################################Clinical data
# load data
icgc.lihc.risk.score <- readRDS(file='/data/icgc_risk_score.rds')
icgc.lihc.cli.data <- readRDS(file='/data/icgc_linc_cli_data.rds')

# data prepare
sig.score <- icgc.lihc.risk.score[, c('icgc_donor_id', 'icgc_specimen_id', 'risk.score', 'risk.categ')]
colnames(sig.score) <- c('icgc_donor_id', 'Tumor_Sample_Barcode', 'risk.score', 'risk.categ')
sig.score$risk.categ <- factor(sig.score$risk.categ, levels=c('low risk', 'high risk'))

icgc.lihc.cli.data$tumor_size_category <- ifelse(icgc.lihc.cli.data$Tumor.size..mm. > 50, '>50', '≤50')
icgc.lihc.cli.data$tumor_size_category <- factor(icgc.lihc.cli.data$tumor_size_category, levels=c('≤50', '>50'))


cli.char <- icgc.lihc.cli.data[, c('icgc_donor_id', 'age_category', 'gender_category', 'Smoking', 
'tumor_size_category', 'virus_infection_category', 'pathologic_stage_category', 'histologic_grade_category', 
'portal_invasion_category', 'hepatic_invasion_category', 'bile_invasion_category', 'fibrosisc_category')]

cli.sig.char <- merge(sig.score, cli.char, by='icgc_donor_id')


# Fisher's Exact Test
p.values <- lapply(colnames(cli.sig.char)[5:15], function(char){
 print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
 print(table(cli.sig.char[, c(char, 'risk.categ')]))
 
 or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
 p.value <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
 return(c(or, p.value))
})

p.values <- as.data.frame(do.call(rbind, p.values))
colnames(p.values) <- c('OR', 'Pvalue')
p.values$Qvalue <- p.adjust(p.values$Pvalue, 'fdr')
rownames(p.values) <- colnames(cli.sig.char)[5:15]



####Plot
library('ggpubr')

# size.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('tumor_size_category', 'risk.categ')])), 
 # "risk.categ", "Freq", fill = "tumor_size_category", color="tumor_size_category", 
 # ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 # title='Tumor size', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

grade.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('histologic_grade_category', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "histologic_grade_category", color="histologic_grade_category", 
 ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='Grade', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

hepatic.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('hepatic_invasion_category', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "hepatic_invasion_category", color="hepatic_invasion_category", 
 ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='Hepatic vein invasion', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

portal.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('portal_invasion_category', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "portal_invasion_category", color="portal_invasion_category", 
 ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='Portal vein invasion', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

ggsave(ggarrange(grade.p, hepatic.p, portal.p, ncol=4, nrow=3, legend = "none"), 
 file='/result/Section4/icgc_cli_risk_com.pdf')



##########################################Mutation data
library('maftools')
lihc.maf <- icgcSimpleMutationToMAF(icgc = '/data/OriginalData/ICGC_LIRI_JP/simple_somatic_mutation.open.tsv.gz', 
 addHugoSymbol = TRUE)
 
lihc.maf$Tumor_Sample_Barcode <- lihc.maf$icgc_specimen_id
lihc.maf <- subset(lihc.maf, Tumor_Sample_Barcode %in% cli.sig.char$Tumor_Sample_Barcode)

mut.maf = read.maf(maf = lihc.maf, clinicalData = cli.sig.char, verbose = TRUE)


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
portal.inva.cols <- cols[7:8]
names(portal.inva.cols) <- c('No', 'Yes')
hepatic.inva.cols <- cols[9:10]
names(hepatic.inva.cols) <- c('No', 'Yes')

 
annocolors = list(risk.categ = risk.cols, portal_invasion_category = portal.inva.cols, 
 pathologic_stage_category = stage.cols, histologic_grade_category = grade.cols, hepatic_invasion_category = hepatic.inva.cols)

pdf('/result/Section4/icgc_mut_oncoplot.pdf')
oncoplot(maf = mut.maf, genes = show_genes, colors = vc_cols, # minMut=9,
  clinicalFeatures = c('risk.categ','pathologic_stage_category', 
  'histologic_grade_category','portal_invasion_category','hepatic_invasion_category'),
  annotationColor = annocolors, sortByAnnotation = TRUE)
dev.off()


# Clinical enrichment analysis
fab.ce = clinicalEnrichment(maf = mut.maf, clinicalFeature = 'risk.categ', minMut = 5)
# fab.ce$pairwise_comparision[fdr < 0.05]

pdf('/result/Section4/icgc_mut_enrich.pdf')
 plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)
dev.off()



# Fisher's Exact Test
lihc.mut.maf <- subset(lihc.maf, !Variant_Classification %in% 
 c("3'Flank", "3'UTR", "5'Flank","5'UTR", "IGR", "Intron", "RNA", "Silent"))

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

        # TP53          RB1       ARID1A         APOB       CTNNB1 
# 6.785544e-03 4.449215e-01 1.000000e+00 3.934958e-01 9.059321e-06 

p.adjust(p.values, method='fdr')
        # TP53          RB1       ARID1A         APOB       CTNNB1 
# 0.0169638599 0.5561518348 1.0000000000 0.5561518348 0.0000452966 





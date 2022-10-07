
lci.risk.score <- readRDS(file='/data/lci_risk_score.rds')

lci.risk.score <- lci.risk.score[, c(3, 5, 	8:18, 24)]
colnames(lci.risk.score) <- c('Affy_GSM',  'Metastasis_risk', 'Gender', 'Age', 'HBV_status', 'ALT', 'Tumor_size', 
 'Multinodular', 'Cirrhosis', 'TNM_staging', 'BCLC_staging', 'CLIP_staging', 'AFP', 'risk.categ')

lci.risk.score$Metastasis_risk <- factor(lci.risk.score$Metastasis_risk, levels=c('low', 'high'))
lci.risk.score$Age <- ifelse(lci.risk.score$Age < 60, '<60', '≥60')
lci.risk.score$ALT <- factor(lci.risk.score$ALT, levels=c('low', 'high'))

lci.risk.score$Tumor_size[lci.risk.score$Tumor_size == '.'] <- NA
lci.risk.score$Tumor_size <- factor(lci.risk.score$Tumor_size, levels=c('small', 'large'))

lci.risk.score$TNM_staging[lci.risk.score$TNM_staging == '.'] <- NA
lci.risk.score$TNM_staging <- ifelse(lci.risk.score$TNM_staging %in% c('I', 'II'), 'I/II', 'III/IV')
lci.risk.score$BCLC_staging[lci.risk.score$BCLC_staging == '.'] <- NA
lci.risk.score$BCLC_staging <- ifelse(lci.risk.score$BCLC_staging %in% c('0', 'A'), '0/A', 'B/C')
lci.risk.score$CLIP_staging[lci.risk.score$CLIP_staging == '.'] <- NA
lci.risk.score$CLIP_staging <- ifelse(lci.risk.score$CLIP_staging %in% c('0', '1'), '0/1', '≥2')
lci.risk.score$AFP[lci.risk.score$AFP == '.'] <- NA

lci.risk.score$AFP <- factor(lci.risk.score$AFP, levels=c('low', 'high'))
lci.risk.score$risk.categ <- factor(lci.risk.score$risk.categ, levels=c('low risk', 'high risk'))


cli.sig.char <- lci.risk.score
# Fisher's Exact Test
p.values <- lapply(colnames(cli.sig.char)[c(2:4, 6:13)], function(char){
 print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
 print(table(cli.sig.char[, c(char, 'risk.categ')]))
 
 or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
 p.value <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
 return(c(or, p.value))
})

p.values <- as.data.frame(do.call(rbind, p.values))
colnames(p.values) <- c('OR', 'Pvalue')
p.values$Qvalue <- p.adjust(p.values$Pvalue, 'fdr')
rownames(p.values) <- colnames(cli.sig.char)[c(2:4, 6:13)]



####Plot
library('ggpubr')

metastasis.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('Metastasis_risk', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "Metastasis_risk", color="Metastasis_risk", ylab='Number of patients', xlab=FALSE, 
 label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9), title='Metastasis risk', 
 lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

tnm.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('TNM_staging', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "TNM_staging", color="TNM_staging", ylab='Number of patients', xlab=FALSE, 
 label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='TNM stage', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

bclc.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('BCLC_staging', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "BCLC_staging", color="BCLC_staging", ylab='Number of patients', xlab=FALSE, 
 label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='BCLC stage', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

clip.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('CLIP_staging', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "CLIP_staging", color="CLIP_staging", ylab='Number of patients', xlab=FALSE, 
 label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='CLIP stage', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

afp.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('AFP', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "AFP", color="AFP", ylab='Number of patients', xlab=FALSE, 
 label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='AFP', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))


ggsave(ggarrange(metastasis.p, tnm.p, bclc.p, clip.p, afp.p, ncol=4, nrow=4, legend = "none"), 
 file='/result/Section4/GSE14520_cli_risk_com.pdf')

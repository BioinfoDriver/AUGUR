
load(file='/data/exp_gene_anno.RData')

gene.exp <- read.csv(file = '/data/OriginalData/Renji_cohort_MR_exp.matrix', 
 sep='\t', row.names = 1, header=TRUE)
colnames(gene.exp) <- gsub('_', '.', colnames(gene.exp))
gene.exp <- gene.exp[, grep('T', colnames(gene.exp))]

rownames(gene.exp) <- gsub('\\.', '-', rownames(gene.exp))

symbol.index <- sapply(rownames(gene.exp), function(symbol){
	index <- match(symbol, exp.gene.anno$Symbol)
	
	if(is.na(index)){
		index <- grep(symbol, exp.gene.anno$Synonyms)[1]	
	}
	
	if(is.na(index)){
		index <- grep(symbol, exp.gene.anno$Symbol_from_nomenclature_authority)[1]	
	} 
	
	return(index)
})

symbol.index <- symbol.index[!is.na(symbol.index) & !duplicated(symbol.index)]
gene.exp <- gene.exp[names(symbol.index), ]
rownames(gene.exp) <- exp.gene.anno$GeneID[symbol.index]


saveRDS(gene.exp, file='/data/renji_norm_gene_exp.rds')

###########################################################
norm.gene.exp <- as.data.frame(t(gene.exp))
norm.gene.exp <- norm.gene.exp %>% tibble::rownames_to_column("SampleID")

# Signature
# signature coefficient
risk.sig <- readRDS(file='/data/curated_lihc_risk_signature.rds')
risk.sig <- risk.sig[, c('GeneID', 'Coefficient', 'PMID')]	
risk.sig$PMID <- paste0('PMID', risk.sig$PMID)
colnames(risk.sig)[1:2] <- c('Gene', 'beta')
risk.sig.list <- split(x=risk.sig, f=risk.sig$PMID)
risk.sig.list <- risk.sig.list[unique(risk.sig$PMID)]


load(file='/data/lasso_cox_res.RData')
active.k.vals <- data.frame(Gene=names(active.k.vals), beta=active.k.vals, PMID='CloSignature')
risk.sig.list <- c(list('CloSignature'=active.k.vals), risk.sig.list)


setwd('/result/Section1/Revise/Renji')
# gene signatures
select.sig <- data.frame(PMID=c('PMID32198063', 'PMID30885723', 'PMID35123387','PMID31335995','PMID31695766',
  'PMID33033585','PMID34676211','PMID34277622','PMID33777771','PMID34975331', 'PMID25666192','PMID22105560', 'CloSignature'),
 SigName=c('Zhang (2020) Genomics', 'Long (2019) EBioMedicine', 'Deng (2022) MolMed', 'Zhang (2019) CancerSci',
  'Long (2019) Theranostics', 'Liu (2020) ComputStructBiotechnolJ',
  'Song (2021) FrontCellDevBiol', 'Xu (2021) FrontCellDevBiol', 'Zhu (2021) FrontOncol',
  'Tang (2022) IntJBiolSci', 'Villa (2016) Gut', 'Kim (2012) Hepatology', 'CloSignature'))


RiskEstiPlot <- function(exp.data, risk.sig, sig.name, cut.off, plot.title, file.name){
  library(forcats)
  library(ggplot2)
  
  # input: title
  # title <- paste0(plot.title, "\nin multiple region LIHC cohort")
  
  # input: signature 
  prognostic_signature <- risk.sig
  
  # input: expression matrix
  gene_matrix <- exp.data
  
  # select signature genes
  inter.sect <- intersect(colnames(gene_matrix), prognostic_signature$Gene)
  
  if(length(prognostic_signature$Gene) > length(inter.sect)){
	print('The following genes are missed:')
	
	print(setdiff(prognostic_signature$Gene, inter.sect))
  }
  
  prognostic_signature <- subset(prognostic_signature, Gene %in% inter.sect)
  tmp <- gene_matrix[, prognostic_signature$Gene] %>% t() %>% as.data.frame()
  colnames(tmp) <- gene_matrix$SampleID
  gene_matrix <- tmp
 
  # calculate RiskScore
  RiskScore <- colSums(gene_matrix * prognostic_signature$beta)
  
  if(sig.name == 'PMID33777771'){
    RiskScore <- exp(RiskScore)
  }
  
  
  RiskScore <- data.frame(SampleID=names(RiskScore), RiskScore=RiskScore)
  RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[1]})
  RiskScore$RegionID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[2]})
  
  
  # define riskscore_thresh
  if(is.numeric(cut.off)){
    riskscore_thresh <- cut.off
  
  }else{
	riskscore_thresh <- median(RiskScore$RiskScore)
  }

  # classify patients as high-risk or low-risk
  RiskScore$RiskScore_bin <- ifelse(RiskScore$RiskScore > riskscore_thresh, "High", "Low")
  
  ## get risk classes
  risk_class <- table(RiskScore$PatientID, RiskScore$RiskScore_bin)
  risk_class <- data.frame(High=as.matrix(risk_class)[,"High"], Low=as.matrix(risk_class)[,"Low"]) %>% tibble::rownames_to_column("PatientID")
  
  # assign as "low", "discordant" or "high"
  risk_class$class <- NA
  risk_class$class <- ifelse(risk_class$High > 0, paste(risk_class$class, "High", sep=""), risk_class$class)
  risk_class$class <- ifelse(risk_class$Low > 0, paste(risk_class$class, "Low", sep=""), risk_class$class)
  risk_class$class <- gsub(x=risk_class$class, pattern="NA", replacement="")
  risk_class$class <- gsub(x=risk_class$class, pattern="HighLow", replacement="Discordant")
  
  # join to RiskScore df
  RiskScore <- dplyr::left_join(x=RiskScore, y=risk_class[,c("PatientID", "class")], by="PatientID")
  RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))
  
  # no. per class
  tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% distinct()
  
  
  # scatter-plot
  sc.plot <- ggplot(RiskScore, aes(x=fct_reorder(PatientID, RiskScore + as.numeric(class), .fun=min), y=RiskScore)) 
  sc.plot <- sc.plot + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(floor(min(RiskScore$RiskScore)), ceiling(max(RiskScore$RiskScore))))
  sc.plot <- sc.plot + theme(axis.text.x = element_text(angle=90, vjust=0.5))
  sc.plot <- sc.plot + geom_hline(yintercept = riskscore_thresh, col="black", lty="dotted")
  sc.plot <- sc.plot + geom_line(col="black")
  sc.plot <- sc.plot + ggtitle(label="", subtitle = paste0(length(unique(RiskScore$PatientID)), " LIVER patients = ", table(tmp$class)["Low"], " low + ", table(tmp$class)["High"], " high + ", table(tmp$class)["Discordant"], " discordant")) + theme(plot.title = element_text(hjust=0.5, face="bold")) + xlab("PatientID") + ylab(plot.title)
  sc.plot <- sc.plot + theme(legend.position = "bottom") + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) #border
  sc.plot <- sc.plot + theme(aspect.ratio = 0.5)
  sc.plot <- sc.plot + geom_point(pch=16, aes(col=class), size=5, alpha=0.5) + scale_color_manual(values = c("#3B4992FF", "azure4", "#EE0000FF")) + theme(legend.position = "none")
  
  
  # calculate percentages for risk classes
  tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
  tmp <- data.frame(table(tmp$class))
  tmp$class <- tmp$Var1 %>% as.character()
  tmp$class[c(1,3)] <- paste0("Concordant\n", tmp$class[c(1,3)], " Risk")
  tmp <- tmp[c(1,3,2),]
  tmp$class <- factor(tmp$class, levels=tmp$class)
  tmp$Perc <- round(tmp$Freq/sum(tmp$Freq)*100)
  
  #  bar-plot
  bar.plot <- ggplot(tmp, aes(x=class, y=Perc, fill=class)) + geom_bar(stat="identity")
  bar.plot <- bar.plot + scale_fill_manual(values = c("#3B4992FF", "#EE0000FF", "gray75"), guide = guide_legend(title=NULL))
  bar.plot <- bar.plot + geom_text(size=5, aes(y = (Perc+2.5), label = paste0(Perc, "%"))) 
  bar.plot <- bar.plot + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(0,round(max(tmp$Perc))+10)) + theme(legend.position = "none", aspect.ratio = 1)
  bar.plot <- bar.plot + xlab("") + ylab("Survival risk classification (%)")
  bar.plot <- bar.plot + ggtitle(label=plot.title)

  pdf(paste0(file.name, '.pdf'))
    print(sc.plot)
    print(bar.plot)
  dev.off()
}


MultiRiskEstiPlot <- function(exp.data, risk.sig.list, cut.off, plot.title, file.name){
  library(forcats)
  library(ggplot2)
  
  # input: title
  # title <- paste0(plot.title, "\nin multiple region LIHC cohort")
  
  discordant.risk <- lapply(names(risk.sig.list), function(sig.name){

    # input: signature 
    prognostic_signature <- risk.sig.list[[sig.name]]
    
    # input: expression matrix
    gene_matrix <- exp.data
    
    # select signature genes
    inter.sect <- intersect(colnames(gene_matrix), prognostic_signature$Gene)
    
    if(length(prognostic_signature$Gene) > length(inter.sect)){
	  print('The following genes are missed:')
	  
	  print(setdiff(prognostic_signature$Gene, inter.sect))
    }
    
    prognostic_signature <- subset(prognostic_signature, Gene %in% inter.sect)
    tmp <- gene_matrix[, prognostic_signature$Gene] %>% t() %>% as.data.frame()
    colnames(tmp) <- gene_matrix$SampleID
    gene_matrix <- tmp
    
    # calculate RiskScore
    RiskScore <- colSums(gene_matrix * prognostic_signature$beta)
    
  	if(sig.name == 'PMID33777771'){
		RiskScore <- exp(RiskScore)
	}	
	
	RiskScore <- data.frame(SampleID=names(RiskScore), RiskScore=RiskScore)
    RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[1]})
    RiskScore$RegionID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[2]})
    
    
    # define riskscore_thresh
    if(is.numeric(cut.off)){
      riskscore_thresh <- cut.off
    
    }else{
	  riskscore_thresh <- median(RiskScore$RiskScore)
    }
    
    # classify patients as high-risk or low-risk
    RiskScore$RiskScore_bin <- ifelse(RiskScore$RiskScore > riskscore_thresh, "High", "Low")
    
    ## get risk classes
    risk_class <- table(RiskScore$PatientID, RiskScore$RiskScore_bin)
    risk_class <- data.frame(High=as.matrix(risk_class)[,"High"], Low=as.matrix(risk_class)[,"Low"]) %>% tibble::rownames_to_column("PatientID")
    
    # assign as "low", "discordant" or "high"
    risk_class$class <- NA
    risk_class$class <- ifelse(risk_class$High > 0, paste(risk_class$class, "High", sep=""), risk_class$class)
    risk_class$class <- ifelse(risk_class$Low > 0, paste(risk_class$class, "Low", sep=""), risk_class$class)
    risk_class$class <- gsub(x=risk_class$class, pattern="NA", replacement="")
    risk_class$class <- gsub(x=risk_class$class, pattern="HighLow", replacement="Discordant")
    
    # join to RiskScore df
    RiskScore <- dplyr::left_join(x=RiskScore, y=risk_class[,c("PatientID", "class")], by="PatientID")
    RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))
    
    
    # calculate percentages for risk classes
    tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
    tmp <- data.frame(table(tmp$class))
    tmp$class <- tmp$Var1 %>% as.character()
    tmp$class[c(1,3)] <- paste0("Concordant ", tmp$class[c(1,3)], " Risk")
    tmp <- tmp[c(1,3,2),]
    tmp$class <- factor(tmp$class, levels=tmp$class)
    # tmp$Perc <- round(tmp$Freq/sum(tmp$Freq)*100)
	tmp$Perc <- as.numeric(sprintf('%0.1f', tmp$Freq/sum(tmp$Freq)*100))
	tmp$SigName <-sig.name 
  
	return(tmp)
  })
  
  discordant.risk <- do.call(rbind, discordant.risk)
  
  #  bar-plot
  # tmp <- subset(discordant.risk, class == 'Discordant')
  # sig.name.order <- tmp$SigName[order(tmp$Perc)]
  
  bar.plot <- ggbarplot(data=discordant.risk, x='SigName', y='Perc', fill='class', color='class', label = TRUE, 
  lab.col = "white", lab.pos = "in", x.text.angle=30, palette = c("#3B4992FF", "#EE0000FF", "gray75"),
  lab.nb.digits=0, xlab = FALSE, ylab="Survival risk classification (%)")#, order=sig.name.order
  
  # ggplot(discordant.risk, aes(x=SigName, y=Perc, fill=class)) + geom_bar(stat="identity", position="fill")
  # bar.plot <- bar.plot + scale_fill_manual(values = c("#3B4992FF", "#EE0000FF", "gray75"), guide = guide_legend(title=NULL))
  # bar.plot <- bar.plot + geom_text(size=5, aes(y = (Perc+2.5), label = paste0(Perc, "%"))) 
  # bar.plot <- bar.plot + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(0,50)) + theme(legend.position = "none", aspect.ratio = 1)
  # bar.plot <- bar.plot + xlab("") + ylab("Survival risk classification (%)")
  # bar.plot <- bar.plot + ggtitle(label=plot.title)

  ggsave(bar.plot, file=paste0(file.name, '.pdf'), width=7)

  return(discordant.risk)
}



for(sig.name in select.sig$PMID){
 
 sig <- risk.sig.list[[sig.name]]
 
 RiskEstiPlot(norm.gene.exp, sig, sig.name, 'median', select.sig$SigName[match(sig.name, select.sig$PMID)], 
 gsub(' ', '_', select.sig$SigName[match(sig.name, select.sig$PMID)]))

}

discor.rate <- MultiRiskEstiPlot(norm.gene.exp, risk.sig.list[select.sig$PMID], 'median', 
 'Tumor sampling bias confounds liver cancer biomarkers', 'select_signature_discordance_rate')


########################### add nault signature
A Hepatocellular Carcinoma 5-Gene Score
five.gene.sig <- cbind(GeneSymbol = c('TAF9','RAMP3','HN1','KRT19','RAN'),
 GeneID = c('6880','10268','51155','3880','5901'))

# Parameter
dqda.centered <- data.frame(GeneSymbol = c('TAF9','RAMP3','HN1','KRT19','RAN'), 
 u.0 = c(-0.16707252, 0.25124668, -0.07128225, 0.4099826, -0.11267109), 
 u.1 = c(0.469731, -0.2855013, 0.4328778, 0.379904, 0.2937122), 
 v.0 = c(0.1501967, 0.2760526, 0.3001276, 1.03919, 0.2766244),
 v.1 = c(0.2989976, 0.2305511, 0.5369335, 0.7997748, 0.4549733))



# Centroid-based (DQDA)
DistanceDQDAPredictor <- function(x, uGood, vGood, uBad, vBad){
 distance.good <- sum(((x-uGood)^2)/vGood) + sum(log(vGood))
 distance.bad <- sum(((x-uBad)^2)/vBad) + sum(log(vBad))
  
 class.label <- ifelse(distance.good < distance.bad, 'Good',
  ifelse(distance.good > distance.bad, 'Bad', 'Undetermined'))

 return(class.label)		 
}

# Stat
DiscordantRisk <- function(risk.class){

 library('dplyr')
 RiskScore <- data.frame(SampleID=names(risk.class), RiskScore_bin=risk.class)
 RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[1]})
 RiskScore$RegionID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[2]})
 
 
 ## get risk classes
 risk_class <- table(RiskScore$PatientID, RiskScore$RiskScore_bin)
 risk_class <- data.frame(High=as.matrix(risk_class)[,"Good"], Low=as.matrix(risk_class)[,"Bad"]) %>% tibble::rownames_to_column("PatientID")
 
 # assign as "low", "discordant" or "high"
 risk_class$class <- NA
 risk_class$class <- ifelse(risk_class$High > 0, paste(risk_class$class, "High", sep=""), risk_class$class)
 risk_class$class <- ifelse(risk_class$Low > 0, paste(risk_class$class, "Low", sep=""), risk_class$class)
 risk_class$class <- gsub(x=risk_class$class, pattern="NA", replacement="")
 risk_class$class <- gsub(x=risk_class$class, pattern="HighLow", replacement="Discordant")
 
 # join to RiskScore df
 RiskScore <- dplyr::left_join(x=RiskScore, y=risk_class[,c("PatientID", "class")], by="PatientID")
 RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))
 
 
 # calculate percentages for risk classes
 tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
 tmp <- data.frame(table(tmp$class))
 tmp$class <- tmp$Var1 %>% as.character()
 tmp$class[c(1,3)] <- paste0("Concordant ", tmp$class[c(1,3)], " Risk")
 tmp <- tmp[c(1,3,2),]
 tmp$class <- factor(tmp$class, levels=tmp$class)
 # tmp$Perc <- round(tmp$Freq/sum(tmp$Freq)*100)
 tmp$Perc <- as.numeric(sprintf('%0.1f', tmp$Freq/sum(tmp$Freq)*100))
 
 return(list(risk_class, tmp))
}


norm.exp.sig <- gene.exp[five.gene.sig[, 2], ]

risk.class.dqda.cen <- apply(t(scale(t(norm.exp.sig), scale = F)), 2, function(y){ # Low=4, High=5, Discordant=2(H3 and H12)
  DistanceDQDAPredictor(y, dqda.centered$u.0, dqda.centered$v.0, dqda.centered$u.1, dqda.centered$v.1)})

DiscordantRisk(risk.class.dqda.cen)

# 1        Low    4  Concordant Low Risk 28.6
# 3       High    5 Concordant High Risk 35.7
# 2 Discordant    5           Discordant 35.7

################################revised
discor.rate <- subset(discor.rate, class == 'Discordant')
discor.rate$sigNum <- sapply(risk.sig.list, nrow)[discor.rate$SigName]


discor.rate <- rbind(discor.rate, c('Discordant', 5, 'Discordant', 35.7, 'PMID23567350', 5))
discor.rate <- dplyr::mutate(discor.rate, Perc=as.numeric(Perc), sigNum = as.numeric(sigNum))


# cor.test(~Perc+sigNum, subset(discor.rate, SigName != 'CloSignature'), method = "spearman") # p-value = 0.2792, rho = 0.1362331
# cor.test(~Perc+sigNum, subset(discor.rate, SigName %in% 
 # setdiff(c(select.sig$PMID, 'PMID23567350'), 'CloSignature')), method = "spearman") # p-value=0.6794, rho=0.1269496


#############################################################

RiskScoreSD <- function(exp.data, risk.sig, sig.name){
    
  # input: signature 
  # input: expression matrix

  
  # select signature genes
  inter.sect <- intersect(rownames(exp.data), names(risk.sig))
  
  if(length(risk.sig) > length(inter.sect)){
	print(sig.name)
	print('The following genes are missed:')
	
	print(setdiff(names(risk.sig), inter.sect))
  }
  
  risk.sig <- risk.sig[inter.sect]
  exp.data <- exp.data[inter.sect, ]

 
  # calculate RiskScore
  RiskScore <- colSums(exp.data * risk.sig)
  
	if(sig.name == 'PMID33777771'){
		RiskScore <- exp(RiskScore)
	}
  
  
  RiskScore <- data.frame(SampleID=names(RiskScore), RiskScore=RiskScore)
  RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[1]})
  
  # SD
  risk.score.sd <- tapply(RiskScore$RiskScore, factor(RiskScore$PatientID), sd)
  return(risk.score.sd)
}

SigScoreSDPlot <- function(risk.sig.list, gene.exp.mat, nault.sig.dist.sd, file.name){
  
  risk.score.sds <- lapply(names(risk.sig.list), function(sig.name){
    risk.score.sd <- RiskScoreSD(gene.exp.mat, risk.sig.list[[sig.name]], sig.name)

    return(risk.score.sd)
  })

  sig.risk.score.sd <- data.frame(risk.score.sd=unlist(risk.score.sds), 
   sig.name = rep(names(risk.sig.list), times=sapply(risk.score.sds, length)))
  
  sig.risk.score.sd <- rbind(sig.risk.score.sd, nault.sig.dist.sd)
  
  label.levels <- names(sort(tapply(sig.risk.score.sd$risk.score.sd, sig.risk.score.sd$sig.name, median)))
  sig.risk.score.sd$sig.name <- factor(sig.risk.score.sd$sig.name, levels=label.levels)
   
  
  print(sort(tapply(sig.risk.score.sd$risk.score.sd, sig.risk.score.sd$sig.name, median)))
  
  # plot 
  library('ggpubr')
  score.sd.plot <- ggboxplot(sig.risk.score.sd, x='sig.name', y='risk.score.sd', fill='sig.name', 
   xlab='Signature', ylab='SD of risk score', legend = "none", add='point', x.text.angle=30)
  
  ggsave(score.sd.plot, file=paste(file.name, '.pdf'), width=7)
}


# Signature
# signature coefficient
load(file='/data/lasso_cox_res.RData')
sig.coef <- readRDS(file='/data/curated_lihc_risk_signature.rds')

risk.sig.list <- lapply(unique(sig.coef$PMID), function(pmid){
  
  risk.sig <- subset(sig.coef, PMID == pmid)$Coefficient
  names(risk.sig) <- subset(sig.coef, PMID == pmid)$GeneID
  
  return(risk.sig)
})
names(risk.sig.list) <- paste0('PMID', unique(sig.coef$PMID))
risk.sig.list <- c(list('CloSignature'=active.k.vals), risk.sig.list)

# RiskScoreSD(losic.gene.exp, risk.sig.list[['Clonal Signature']])

risk.sig <- c('CloSignature', 'PMID32198063', 'PMID30885723', 'PMID35123387', 'PMID31335995', 'PMID31695766',
 'PMID33033585', 'PMID34676211', 'PMID34277622', 'PMID33777771', 'PMID34975331', 'PMID25666192','PMID22105560')
risk.sig.list <- risk.sig.list[risk.sig] 

########################### add nault signature
# A Hepatocellular Carcinoma 5-Gene Score
A Hepatocellular Carcinoma 5-Gene Score
five.gene.sig <- cbind(GeneSymbol = c('TAF9','RAMP3','HN1','KRT19','RAN'),
 GeneID = c('6880','10268','51155','3880','5901'))

# Parameter
dqda.centered <- data.frame(GeneSymbol = c('TAF9','RAMP3','HN1','KRT19','RAN'), 
 u.0 = c(-0.16707252, 0.25124668, -0.07128225, 0.4099826, -0.11267109), 
 u.1 = c(0.469731, -0.2855013, 0.4328778, 0.379904, 0.2937122), 
 v.0 = c(0.1501967, 0.2760526, 0.3001276, 1.03919, 0.2766244),
 v.1 = c(0.2989976, 0.2305511, 0.5369335, 0.7997748, 0.4549733))


# Centroid-based (DQDA)
DistanceDQDAPredictorScore <- function(x, uGood, vGood, uBad, vBad){
 distance.good <- sum(((x-uGood)^2)/vGood) + sum(log(vGood))
 distance.bad <- sum(((x-uBad)^2)/vBad) + sum(log(vBad))
  
 class.score <- c(distance.good, distance.bad)
 names(class.score) <- c('Good', 'Bad')

 return(class.score)		 
}


norm.exp.sig <- gene.exp[five.gene.sig[, 2], ]

renji.dist <- apply(t(scale(t(norm.exp.sig), scale = F)), 2, function(y){
  DistanceDQDAPredictorScore(y, dqda.centered$u.0, dqda.centered$v.0, dqda.centered$u.1, dqda.centered$v.1)})

renji.dist <- as.data.frame(t(renji.dist))
renji.dist$PatientID <- sapply(X=rownames(renji.dist), FUN=function(x) {unlist(strsplit(as.character(x), split="\\."))[1]})

renji.dist.sd <- data.frame(risk.score.sd=tapply(renji.dist$Good, renji.dist$PatientID, sd), sig.name = 'PMID23567350')

########################### add nault signature

setwd('/result/Section1/Revise/Renji')
SigScoreSDPlot(risk.sig.list, gene.exp, renji.dist.sd, 'sig_score_sd_in_renji')

# PMID31695766 CloSignature PMID30885723 PMID34975331 PMID34277622 PMID35123387 PMID34676211 PMID33033585 PMID31335995 PMID33777771 
   # 0.0897240    0.1706793    0.1824075    0.1872940    0.2052390    0.2214055    0.2389983    0.3154026    0.6594846    0.7472416 
# PMID32198063 PMID25666192 PMID22105560 PMID23567350 
   # 0.8268226    1.0015902    2.4488576    5.3372133


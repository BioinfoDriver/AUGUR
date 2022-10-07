
# load data
# Gene expression from Losic et al.
losic.gene.exp <- readRDS(file='/data/losic_vst_norm_geneExp.rds')

losic.gene.exp <- losic.gene.exp[, c(
"H1.a", "H1.b", 
"H2.a1", "H2.b", "H2.c", "H2.d", "H2.e", #  "H2.a2"
"H3.a", "H3.b", 
"H4.a", "H4.b", "H4.c", "H4.d", "H4.e", 
"H6.a", "H6.b", 
"H7.a", "H7.b", "H7.c", "H7.d", "H7.e", 
"H8.a", "H8.b", "H8.c", 
"H9.a", "H9.b", "H9.c", "H9.d", "H9.e", "H9.f", 
"H10.a", "H10.b", "H10.c", "H10.d", "H10.e", 
"H11.a", "H11.b", 
"H12.a", "H12.b", "H12.c", "H12.d", "H12.e")]
colnames(losic.gene.exp)[3] <- "H2.a"

losic.gene.exp <- as.data.frame(t(losic.gene.exp))
losic.gene.exp <- losic.gene.exp %>% tibble::rownames_to_column("SampleID")


# Gene expression from Shi et al.
shi.gene.exp <- readRDS(file='/data/shi.lihc.gene.exp.rds')
shi.gene.exp <- shi.gene.exp[, c(
"H1.1", "H1.2", "H1.3", "H1.4", "H1.5",
"H2.1", "H2.2", "H2.3", "H2.4", "H2.5", 
"H3.1", "H3.2", "H3.3", "H3.4", "H3.5", 
"H4.1", "H4.2", "H4.3", "H4.4", "H4.5", 
"H5.1", "H5.2", "H5.3", "H5.4", "H5.5")]

shi.gene.exp <- as.data.frame(t(shi.gene.exp))
shi.gene.exp <- shi.gene.exp %>% tibble::rownames_to_column("SampleID")


# plot
# title: e.g. 'Huang (2020) Front Oncol signature'
RiskEstiPlot <- function(exp.data, risk.sig, cut.off, plot.title, file.name){
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


setwd('/result/Section1')
# gene signatures
select.sig <- data.frame(PMID=c('PMID32198063', 'PMID30885723', 'PMID35123387','PMID31335995','PMID31695766',
  'PMID33033585','PMID34676211','PMID34277622','PMID33777771','PMID34975331', 'PMID25666192','PMID22105560'),
 SigName=c('Zhang (2020) Genomics', 'Long (2019) EBioMedicine', 'Deng (2022) MolMed', 'Zhang (2019) CancerSci',
  'Long (2019) Theranostics', 'Liu (2020) ComputStructBiotechnolJ',
  'Song (2021) FrontCellDevBiol', 'Xu (2021) FrontCellDevBiol', 'Zhu (2021) FrontOncol',
  'Tang (2022) IntJBiolSci', 'Villa (2016) Gut', 'Kim (2012) Hepatology'))


for(index in seq(nrow(select.sig))){
 
 sig <- risk.sig.list[[select.sig$PMID[index]]]
 
 RiskEstiPlot(losic.gene.exp, sig, 'median', select.sig$SigName[index], 
 gsub(' ', '_', select.sig$SigName[index]))

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

discor.rate <- MultiRiskEstiPlot(losic.gene.exp, risk.sig.list[select.sig$PMID], 'median', 
 'Tumor sampling bias confounds liver cancer biomarkers', 'select_signature_discordance_rate')


discor.rate[seq(from = 1, to = 65*3, by = 3), ]
discor.rate[seq(from = 2, to = 65*3, by = 3), ]
discor.rate[seq(from = 3, to = 65*3, by = 3), ]


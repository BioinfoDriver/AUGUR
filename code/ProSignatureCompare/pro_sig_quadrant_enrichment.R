
# load data
# Gene expression from Losic et al.
losic.gene.exp <- readRDS(file='/data/losic_vst_norm_geneExp.rds')
load(file='/data/exp_gene_anno.RData')

losic.gene.exp <- losic.gene.exp[as.character(subset(exp.gene.anno, !is.na(LosicExpGene))$GeneID), ]

losic.pats <- list(
H2 = c("H2.a1", "H2.b", "H2.c", "H2.d", "H2.e"), 
H4 = c("H4.a", "H4.b", "H4.c", "H4.d", "H4.e"), 
H7 = c("H7.a", "H7.b", "H7.c", "H7.d", "H7.e"), 
H9 = c("H9.a", "H9.b", "H9.c", "H9.d", "H9.e", "H9.f"), 
H10 = c("H10.a", "H10.b", "H10.c", "H10.d", "H10.e"),
H12 = c("H12.a", "H12.b", "H12.c", "H12.d", "H12.e"))



IntraHetScore <- function(exp.dat, sample.list, method){
	
	ith.scores <- lapply(sample.list, function(sample.set){
		pat.exp.dat <- exp.dat[, sample.set]
		
		if(method == 'sd'){
		  ith.score <- apply(pat.exp.dat, 1, sd)
		}else if(method == 'mad'){
		  ith.score <- apply(pat.exp.dat, 1, mad)
		}else{
		  ith.score <- apply(pat.exp.dat, 1, function(x) sd(x)/mean(x))
		}
		
		return(ith.score)
	})
	ith.scores <- do.call(cbind, ith.scores)
	return(rowMeans(ith.scores))
}


InterHetScore <- function(exp.dat, sample.list, method){
	
	set.seed(41)
	sample.list <- lapply(1:10, function(i){
		sample.set <- sapply(sample.list, function(sample.set) sample(sample.set, 1))
		return(sample.set)
	})
		
	ith.scores <- lapply(sample.list, function(sample.set){
		pat.exp.dat <- exp.dat[, sample.set]
		
		if(method == 'sd'){
		  ith.score <- apply(pat.exp.dat, 1, sd)
		}else if(method == 'mad'){
		  ith.score <- apply(pat.exp.dat, 1, mad)
		}else{
		  ith.score <- apply(pat.exp.dat, 1, function(x) sd(x)/mean(x))
		}
		
		return(ith.score)
	})
	ith.scores <- do.call(cbind, ith.scores)
	return(rowMeans(ith.scores))
}



losic.ith.score <- data.frame(losic.intra.sd = IntraHetScore(losic.gene.exp, losic.pats, 'sd'),
 losic.inter.sd = InterHetScore(losic.gene.exp, losic.pats, 'sd'))

# load(file='/data/intra.inter.ith.score.RData')

within.ith <- ifelse(losic.ith.score$losic.intra.sd > mean(losic.ith.score$losic.intra.sd), 'top', 'bottom')
between.ith <- ifelse(losic.ith.score$losic.inter.sd > mean(losic.ith.score$losic.inter.sd), 'right', 'left')
losic.ith.score$quadrant <- paste(within.ith, between.ith, sep='_')


# Signature
curated.sig <- readRDS(file='/data/curated_lihc_risk_signature.rds')
evaluated.sig <- c('32198063', '30885723', '35123387', '31335995', '31695766',
 '33033585', '34676211', '34277622', '33777771', '34975331', '25666192','22105560')
evaluated.sig <- subset(curated.sig, PMID %in% evaluated.sig)


evaluated.sig$quadrant <- losic.ith.score[evaluated.sig$GeneID, ]$quadrant


# plot
sig.qua.perc <- as.data.frame(prop.table(table(evaluated.sig[, c('PMID', 'quadrant')]), margin=1))
sig.qua.perc$quadrant <- factor(sig.qua.perc$quadrant,levels=c("top_left", "bottom_left", "top_right", "bottom_right"))
sig.qua.perc$PMID <- factor(sig.qua.perc$PMID,levels=c('34277622', '31695766', '34975331', '33777771', '30885723',
 '35123387', '33033585', '34676211', '31335995', '32198063', '25666192', '22105560'))

library('ggpubr')
sig.qua.perc.plot <- ggbarplot(sig.qua.perc, 'PMID', 'Freq', fill = "quadrant", color = "quadrant", x.text.angle=15, 
 palette=c("firebrick1", "darkorchid2", "gold1", "turquoise2"), xlab='Signature', ylab='Percentage of gene')

ggsave(sig.qua.perc.plot, filename='/result/Section5/sig_qua_perc.pdf')


# plot
sig.stat <- subset(evaluated.sig, !duplicated(GeneID) & !is.na(quadrant))
# sig.stat <- subset(evaluated.sig, !duplicated(GeneID) & !is.na(quadrant) & PMID != '22105560')

sig.stat <- prop.table(cbind(table(sig.stat$quadrant), 
 round(nrow(sig.stat)*table(losic.ith.score$quadrant)/sum(nrow(losic.ith.score)))), margin=2)

sig.stat <- as.data.frame(sig.stat)
colnames(sig.stat) <- c('Observed', 'Expeted')

sig.stat <- tibble::rownames_to_column(sig.stat, var = "quadrant")
sig.stat <- reshape2::melt(sig.stat, id='quadrant')
sig.stat$quadrant <- factor(sig.stat$quadrant,levels=c("top_left", "bottom_left", "top_right", "bottom_right"))

sig.qua.stat.plot <- ggbarplot(sig.stat, 'variable', 'value', fill = "quadrant", color = "quadrant", 
 x.text.angle=15, palette=c("firebrick1", "darkorchid2", "gold1", "turquoise2"), xlab=FALSE, 
 ylab='Percentage of gene', label=TRUE, lab.pos='in')

ggsave(sig.qua.stat.plot, filename='/result/Section5/sig_qua_stat.pdf')

# focal on l4321 genes
# P value
# left vs right
chisq.test(matrix(c(45, 85, 81, 49), 2)) # p=1.404e-05
# chisq.test(matrix(c(29, 37, 42, 25), 2)) # delete 22105560, p=0.04625

# top vs bottom
chisq.test(matrix(c(100, 30, 49, 81), 2)) # p=3.632e-10
# chisq.test(matrix(c(44, 22, 25, 42), 2)) # delete 22105560, p=0.001309

# topleft
chisq.test(matrix(c(21, 109, 13, 117), 2)) # p=0.1979
# chisq.test(matrix(c(11, 55, 7, 60), 2)) # delete 22105560, p=0.4268

# topright
chisq.test(matrix(c(79, 51, 36, 94), 2)) # p=1.567e-07
# chisq.test(matrix(c(33, 33, 18, 49), 2)) # delete 22105560, p=0.01031

# bottomleft
chisq.test(matrix(c(24, 106, 68, 62), 2)) # p=2.446e-08
# chisq.test(matrix(c(18, 48, 35, 32), 2)) # delete 22105560, p=0.005723

# bottomright
chisq.test(matrix(c(6, 124, 13, 117), 2)) # p=0.1528
# chisq.test(matrix(c(4, 62, 7, 60), 2)) # delete 22105560, p=0.5461



# focal on 19170 genes
# P value
# left vs right
chisq.test(matrix(c(31, 101, 79, 54), 2)) # p=6.331e-09
# chisq.test(matrix(c(19, 48, 40, 27), 2)) # delete 22105560, p=0.0005007

# top vs bottom
chisq.test(matrix(c(114, 18, 53, 80), 2)) # p=1.209e-14
# chisq.test(matrix(c(55, 12, 27, 40), 2)) # delete 22105560, p=1.698e-06

# topleft
chisq.test(matrix(c(17, 115, 13, 120), 2)) # p=0.5461
# chisq.test(matrix(c(9, 58, 7, 60), 2)) # delete 22105560, p=0.7899

# topright
chisq.test(matrix(c(97, 35, 40, 93), 2)) # p=3.713e-12
# chisq.test(matrix(c(46, 21, 20, 47), 2)) # delete 22105560, p=1.562e-05

# bottomleft
chisq.test(matrix(c(14, 118, 66, 67), 2)) # p=1.169e-11
# chisq.test(matrix(c(10, 57, 33, 37), 2)) # delete 22105560, p=0.0001054

# bottomright
chisq.test(matrix(c(4, 128, 14, 119), 2)) # p=0.0292
# chisq.test(matrix(c(2, 65, 7, 60), 2)) # delete 22105560, p=0.1674

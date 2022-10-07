
# load data
# Gene expression from Losic et al.
losic.gene.exp <- readRDS(file='/data/losic_vst_norm_geneExp.rds')
load(file='/data/exp_gene_anno.RData')

losic.gene.exp <- losic.gene.exp[as.character(subset(exp.gene.anno, !is.na(LosicExpGeneFilter))$GeneID), ]

losic.pats <- list(
# H1 = c("H1.a", "H1.b"),
H2 = c("H2.a1", "H2.b", "H2.c", "H2.d", "H2.e"), # "H2.a2", 
# H3 = c("H3.a", "H3.b"),
H4 = c("H4.a", "H4.b", "H4.c", "H4.d", "H4.e"), 
# H6 = c("H6.a", "H6.b"),
H7 = c("H7.a", "H7.b", "H7.c", "H7.d", "H7.e"), 
# H8 = c("H8.a", "H8.b", "H8.c"),
H9 = c("H9.a", "H9.b", "H9.c", "H9.d", "H9.e", "H9.f"), 
H10 = c("H10.a", "H10.b", "H10.c", "H10.d", "H10.e"),
# H11 = c("H11.a", "H11.b"),
H12 = c("H12.a", "H12.b", "H12.c", "H12.d", "H12.e"))


# Gene expression from Shi et al.
shi.gene.exp <- readRDS(file='/data/shi.lihc.gene.exp.rds')

shi.pats <- list(
H1 = c("H1.1", "H1.2", "H1.3", "H1.4", "H1.5"),
H2 = c("H2.1", "H2.2", "H2.3", "H2.4", "H2.5"), 
H3 = c("H3.1", "H3.2", "H3.3", "H3.4", "H3.5"), 
H4 = c("H4.1", "H4.2", "H4.3", "H4.4", "H4.5"), 
H5 = c("H5.1", "H5.2", "H5.3", "H5.4", "H5.5"))


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
losic.intra.mad = IntraHetScore(losic.gene.exp, losic.pats, 'mad'),
losic.intra.cv = IntraHetScore(losic.gene.exp, losic.pats, 'cv'),
losic.inter.sd = InterHetScore(losic.gene.exp, losic.pats, 'sd'),
losic.inter.mad = InterHetScore(losic.gene.exp, losic.pats, 'mad'),
losic.inter.cv = InterHetScore(losic.gene.exp, losic.pats, 'cv'))


shi.ith.score <- data.frame(shi.intra.sd = IntraHetScore(shi.gene.exp, shi.pats, 'sd'),
shi.intra.mad = IntraHetScore(shi.gene.exp, shi.pats, 'mad'),
shi.intra.cv = IntraHetScore(shi.gene.exp, shi.pats, 'cv'),
shi.inter.sd = InterHetScore(shi.gene.exp, shi.pats, 'sd'),
shi.inter.mad = InterHetScore(shi.gene.exp, shi.pats, 'mad'),
shi.inter.cv = InterHetScore(shi.gene.exp, shi.pats, 'cv'))


ith.score <- merge(losic.ith.score, shi.ith.score, by='row.names')
ith.score <- tibble::column_to_rownames(ith.score, var = "Row.names")

save(losic.ith.score, shi.ith.score, ith.score, file='/data/intra.inter.ith.score.RData')

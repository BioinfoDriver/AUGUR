
# load data
univ.cox.pvalue <- readRDS(file='/data/tcga.lihc.univ.cox.pvalue.rds')

load(file='/data/intra.inter.ith.score.RData')
load(file='/data/icgc_lihc_diff_exp.RData')

# differential gene expression selecction
diff.exp.gene <- subset(icgc.paired.diff.res, abs(log2FoldChange) > 1 & padj < 0.05) # 2024(up), 1221(down)
diff.exp.gene.list <- rownames(diff.exp.gene) # 3245


# clonal gene expression selection
within.ith <- ifelse(losic.ith.score$losic.intra.sd > mean(losic.ith.score$losic.intra.sd), 'top', 'bottom')
between.ith <- ifelse(losic.ith.score$losic.inter.sd > mean(losic.ith.score$losic.inter.sd), 'right', 'left')
losic.ith.score$quadrant <- paste(within.ith, between.ith, sep='_')

clonal.exp.gene.list <- rownames(subset(losic.ith.score, quadrant == 'bottom_right')) # 1477


# prognostic gene selecction
# univ.cox.pvalue <- p.adjust(univ.cox.pvalue, method='fdr')
prog.gene.list <- names(univ.cox.pvalue[univ.cox.pvalue < 0.05]) # 3495
# prog.gene.list <- names(univ.cox.pvalue[univ.cox.pvalue < 0.1]) # 2432

# candidate prognostic signature
candi.prog.sig.gene <- Reduce(intersect, list(diff.exp.gene.list, clonal.exp.gene.list, prog.gene.list))

saveRDS(candi.prog.sig.gene, file='/data/candi_prog_sig_gene.rds')


######Plot
library("VennDiagram")
library("RColorBrewer")

myCol <- brewer.pal(3, "Pastel2")
p <- venn.diagram(x = list(diff.exp.gene.list, clonal.exp.gene.list, prog.gene.list),
	category.names = c("DE Set" , "Q4 Set" , "Surv Set"), filename=NULL, output=TRUE, 
	fill=myCol, fontface = "bold", fontfamily = "sans")

pdf("/result/Section2/candidate_gene_venn_diagramm.pdf")
 grid.draw(p)
dev.off()


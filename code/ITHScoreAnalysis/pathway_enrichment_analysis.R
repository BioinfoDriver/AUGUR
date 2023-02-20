
# load data
load(file='/data/intra.inter.ith.score.RData')

# RNA intra- and intertumor heterogeneity
ii.ith.score <- losic.ith.score[, c('losic.intra.sd', 'losic.inter.sd')]
colnames(ii.ith.score) <- c('intra.score', 'inter.score')
within.ith <- ifelse(ii.ith.score$intra.score > mean(ii.ith.score$intra.score), 'top', 'bottom')
between.ith <- ifelse(ii.ith.score$inter.score > mean(ii.ith.score$inter.score), 'right', 'left')
ii.ith.score$quadrant <- paste(within.ith, between.ith, sep='_')


background.genes <- rownames(ii.ith.score)
q1.genes <- rownames(ii.ith.score)[ii.ith.score$quadrant == 'top_right'] # Q1
q2.genes <- rownames(ii.ith.score)[ii.ith.score$quadrant == 'top_left'] # Q2
q3.genes <- rownames(ii.ith.score)[ii.ith.score$quadrant == 'bottom_left'] # Q3
q4.genes <- rownames(ii.ith.score)[ii.ith.score$quadrant == 'bottom_right'] # Q4


# Pathway enrichment analysis
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)

# GO over-representation analysis
q1.enrich.go <- enrichGO(
  gene = q1.genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  # universe,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE)

q2.enrich.go <- enrichGO(
  gene = q2.genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  # universe,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE)
  
q3.enrich.go <- enrichGO(
  gene = q3.genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  # universe,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE)
  
q4.enrich.go <- enrichGO(
  gene = q4.genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  # universe,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE)

write.table(q1.enrich.go, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE, file='/data/Q1.enrich.go.txt')

write.table(q4.enrich.go, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE, file='/data/Q4.enrich.go.txt')

 
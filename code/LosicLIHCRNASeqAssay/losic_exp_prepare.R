

# Importing transcript abundance with tximport
library(DESeq2)
library(tximport)
library(readr)

setwd('/data/OriginalData/Losic_Nat_Commun_2020_LiverMultiregionExp')
sam.info <- read.csv(file='E-MTAB-5905.sdrf.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

ena.runs <- sam.info$Comment.ENA_RUN.
tran.files <- file.path('/data/OriginalData/RSEM_exp/res_refseq2', ena.runs, paste0(ena.runs, '.isoforms.results'))
names(tran.files) <- sam.info$Extract.Name

tx2gene <- read.csv(file=tran.files[1], header=TRUE, sep='\t', stringsAsFactors=FALSE)[, 1:2]
tran.exp <- tximport(files = tran.files, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)


save(tran.exp, file='/data/losic.refseq.gene.exp.RData')

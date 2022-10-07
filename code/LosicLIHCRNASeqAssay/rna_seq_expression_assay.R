
# wget -i filelist.txt
setwd('/data/OriginalData/Losic_Nat_Commun_2020_LiverMultiregionExp')
sam.info <- read.csv(file='E-MTAB-5905.sdrf.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

# Run RSEM
RsemExpCal <- function(ena.run){
	library('R.utils')
	print(ena.run)
	
	decompressFile(filename=paste0('./RNA-seq/', ena.run, '.fastq.gz'), 
	 destname=paste0('./fastq/', ena.run, '.fastq'), FUN=gzfile, ext='gz', remove=FALSE)

	# dir.create(file.path('./res', ena.run))

	# rsem-calculate-expression --phred33-quals -p 8 --append-names --no-bam-output \
		# ./fastq/ERR2039934.fastq \
		# ./ref/human_ensembl \
		# ./res/ERR2039934/ERR2039934

	# command <- paste('rsem-calculate-expression --phred33-quals -p 8 --append-names --no-bam-output',
	 # paste('./fastq/', ena.run, '.fastq', sep=''), './ref/human_ensembl', file.path('./res', ena.run, ena.run))
	
	dir.create(file.path('./res_refseq2', ena.run))
	
	command <- paste('rsem-calculate-expression --phred33-quals -p 8 --append-names --no-bam-output',
	 paste('./fastq/', ena.run, '.fastq', sep=''), './refseq_ref2/human_refseq', file.path('./res_refseq2', ena.run, ena.run))
	
	
	system(command)

	file.remove(paste0('./fastq/', ena.run, '.fastq'))
	return(NULL)
}

setwd('/data/OriginalData/RSEM_exp')
ena.runs <- sam.info$Comment.ENA_RUN.[1:15]
for(ena.run in ena.runs) try(RsemExpCal(ena.run))

ena.runs <- sam.info$Comment.ENA_RUN.[16:30]
for(ena.run in ena.runs) try(RsemExpCal(ena.run))

ena.runs <- sam.info$Comment.ENA_RUN.[31:45]
for(ena.run in ena.runs) try(RsemExpCal(ena.run))

ena.runs <- sam.info$Comment.ENA_RUN.[46:60]
for(ena.run in ena.runs) try(RsemExpCal(ena.run))

ena.runs <- sam.info$Comment.ENA_RUN.[61:62]
for(ena.run in ena.runs) try(RsemExpCal(ena.run))



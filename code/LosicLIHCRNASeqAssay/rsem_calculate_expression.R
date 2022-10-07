############################################# RSEM
# http://deweylab.github.io/RSEM/README.html

################# Usage
#### I. Preparing Reference Sequences ####
# RSEM can extract reference transcripts from a genome if you provide it with gene annotations in a GTF/GFF3 file.
# To prepare the reference sequences, you should run the rsem-prepare-reference program.
# http://deweylab.github.io/RSEM/rsem-prepare-reference.html

# Build RSEM references using RefSeq, Ensembl, or GENCODE annotations
# RefSeq and Ensembl are two frequently used annotations. For human and mouse, GENCODE annotaions are also available. 
# In this section, we show how to build RSEM references using Ensembl annotation. 
# Note that it is important to pair the genome with the annotation file for each annotation source. 
# In addition, we recommend users to use the primary assemblies of genomes. 

# Download and decompress the human genome and GTF files to your working directory:
# e.g. /pub5/xiaoyun/Jobs/J22/RSEM_exp/genome_gtf
# http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz


# Then use the following command to build RSEM references:
# rsem-prepare-reference [options] reference_fasta_file(s) reference_name

# @ reference_fasta_file(s)
	# Either a comma-separated list of Multi-FASTA formatted files OR a directory name. 
	# If a directory name is specified, RSEM will read all files with suffix ".fa" or ".fasta" in this directory. 

# @ reference name
	# The name of the reference used. 
	# RSEM will generate several reference-related files that are prefixed by this name. 
	# This name can contain path information (e.g. '/ref/mm9').

# @ --gtf <file>
	# If this option is on, RSEM assumes that 'reference_fasta_file(s)' contains the sequence of a genome, 
	# and will extract transcript reference sequences using the gene annotations specified in <file>, which should be in GTF format.

# e.g.
# cd /data/OriginalData/RSEM_exp/
# rsem-prepare-reference --gtf ./genome_gtf/Homo_sapiens.GRCh38.104.gtf \
               # --bowtie \
               # ./genome_gtf/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
               # ./ref/human_ensembl
# OR
# rsem-prepare-reference --gtf ./refseq_gtf/GCF_000001405.25_GRCh37.p13_genomic.gtf \
               # --bowtie \
			   # ./refseq_gtf/GCF_000001405.25_GRCh37.p13_genomic.primary_assembly.fna \
			   # ./refseq_ref/human_refseq
# OR
# rsem-prepare-reference --gff3 ./refseq_gtf/GCF_000001405.25_GRCh37.p13_genomic.gff \
               # --trusted-sources BestRefSeq,Curated\ Genomic \
			   # --gff3-RNA-patterns mRNA \
			   # --bowtie \
			   # ./refseq_gtf/GCF_000001405.25_GRCh37.p13_genomic.primary_assembly.fna \
			   # ./refseq_ref2/human_refseq


# *Default: the path to RSEM executables is assumed to be in the user's PATH environment variable
# *Default: the path to Bowtie executables is assumed to be in the user's PATH environment variable



#### II. Calculating Expression Values ####
# To calculate expression values, you should run the rsem-calculate-expression program.
# http://deweylab.github.io/RSEM/rsem-calculate-expression.html

# rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name 
# rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name 

# @ upstream_read_files(s)
#	Comma-separated list of files containing single-end reads or upstream reads for paired-end data.

# @ downstream_read_file(s)
#	Comma-separated list of files containing downstream reads which are paired with the upstream reads. 

# @ reference_name
# 	The name of the reference used. The user must have run 'rsem-prepare-reference' with this reference_name before running this program.

# @ sample_name
# 	The name of the sample analyzed. All output files are prefixed by this name (e.g., sample_name.genes.results)

# In its default mode, this program aligns input reads against a reference transcriptome with Bowtie and calculates expression values using the alignments. 
# RSEM assumes the data are single-end reads with quality scores, unless the '--paired-end' or '--no-qualities' options are specified.

# e.g.
# cd /data/OriginalData/RSEM_exp/
# rsem-calculate-expression --phred33-quals -p 8 --append-names --no-bam-output \
	# ./fastq/ERR2039934.fastq \
	# ./ref/human_ensembl \
    # ./res/ERR2039934/ERR2039934
	
# @ -p/--num-threads <int>
	# Number of threads to use. Both Bowtie/Bowtie2, expression estimation and 'samtools sort' will use this many threads.

# @ --phred33-quals
	# Input quality scores are encoded as Phred+33. This option is used by Bowtie, Bowtie 2 and HISAT2. (Default: on)

# @ --append-names
	# If gene_name/transcript_name is available, append it to the end of gene_id/transcript_id (separated by '_') in files 
	# 'sample_name.isoforms.results' and 'sample_name.genes.results'.

# @ --no-bam-output
	# Do not output any BAM file.

#### NOTES ####
# For single-end data, it is strongly recommended that the user provide the fragment length distribution parameters (--fragment-length-mean and --fragment-length-sd). 
# For paired-end data, RSEM will automatically learn a fragment length distribution from the data.

# @	--fragment-length-mean <double>
	# (single-end data only) The mean of the fragment length distribution, which is assumed to be a Gaussian.

# @ --fragment-length-sd <double>
	# (single-end data only) The standard deviation of the fragment length distribution, which is assumed to be a Gaussian.

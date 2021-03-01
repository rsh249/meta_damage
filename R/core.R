# R script to align sequence reads to reference database and identify ancient DNA damage patterns 

# set up path to data 
# set up path to reference
source('R/paths.R')

# other ways of downloading packrat midden data
# on the command line/terminal
# fastq-dump example: 
# fastq-dump SRR7774472
# system('fastq-dump SRR7774472')

# is there an R solution using NCBI Entrez tools or sratoolkit?


# getting the reference data
# ftp://ftp.ncbi.nlm.nih.gov/
download.file('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz',
              'mito_reference.fna.gz')
ref = 'mito_reference.fna.gz' # might have to unzip with gunzip


# Consider overlapping reads and filtering for replicates before aligning with bwa

#Quality control: low quality, adapter trimming, duplicates

# align with bwa
index_bwa = paste('bwa index', ref)
system(index_bwa) # index reference
run_bwa = paste('bwa mem -Y -I 0 -L 1024 -E 7 -t 32', ref, data, '> align.sam')
system(run_bwa) # run alignment


# samtools index and convert

run_samtools = paste(
'samtools view -S -b align.sam > align.bam \n
samtools view -b -F4 align.bam > align.final.bam \n
samtools sort align.final.bam > align.sort.bam \n
samtools index align.sort.bam \n')
system(run_samtools)

#### Set paths

system('git clone https://github.com/rsh249/PMDtools') #python3 compatible code?

#pmdtools deamination mapping
# run_pmdtools = paste('samtools calmd -b align.sort.bam', ref, '| samtools view -h - | python3 ../PMDtools/pmdtools.0.60.py --threshold 0 --header --deamination --range 125 --CpG > out.CpG')
# system(run_pmdtools)
run_pmdtools2 = paste('samtools view -h align.sort.bam | python3 ../PMDtools/pmdtools.0.60.py --deamination --range 125 > out2.deam') #as in Moore etal. 2020
system(run_pmdtools2)


#### FILTERING ######
#pmdtools filtering on PMD score 0
run_pmdtools_filter = paste('samtools view -h align.sort.bam | python3 ../PMDtools/pmdtools.0.60.py --threshold 0 --header | samtools view -Sb - > filter.bam') 
system(run_pmdtools_filter)

# bam to fastq
run_bam2fastq = paste('samtools fastq -@ 8 filter.bam -1 filter.R1.fastq -2 filter.R2.fastq')
system(run_bam2fastq)

# plot
library(ShortRead)
library(ggplot2)
filter_fastq = readFastq('.', pattern='*.fastq')
filter_reads = sread(filter_fastq)
filter_readlengths = filter_reads@ranges@width
qscores = quality(filter_fastq)
numqscores = as(qscores, "matrix") # converts to numeric scores automatically
avgscores = apply(numqscores, 1, mean, na.rm=TRUE)
avgscores = as.data.frame(avgscores)
ggplot(avgscores) +
  geom_histogram(aes(x=avgscores), binwidth=0.2) +
  theme_linedraw() +
  xlab('Quality Score') +
  ggtitle('Per Read Average Quality')


# overlap reads and plot length distribution


#out2.deam

# R Script to Align Sequence Reads to Reference Database and Identify Ancient DNA Damage Patterns 

system('fastq-dump SRR7774472')    #### is there an R solution using NCBI Entrez tools or sratoolkit?

# set up path to data 
# set up path to reference
source('R/paths.R')

# other ways of downloading packrat midden data
# on the command line/terminal
# fastq-dump example: 
# fastq-dump SRR7774472
# system('fastq-dump SRR7774472')

# is there an R solution using NCBI Entrez tools or sratoolkit?

# Getting the Reference Data
download.file('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz', 'plastid.1.1.genomic.fna.gz')
download.file('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz', 'plastid.2.1.genomic.fna.gz')
download.file('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz', 'plastid.3.1.genomic.fna.gz')
download.file('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz', 'plastid.4.1.genomic.fna.gz')
ref = 'plastid_all.fna'
data='SRR7774472.fastq'

# Alignment with BWA
index_bwa = paste('bwa index', ref)
system(index_bwa) # index reference

run_bwa = paste('bwa mem -Y -I 0 -L 1024 -E 7 -t 32', ref, data, '> align.sam')
system(run_bwa) # run alignment

#Samtools Index and Convert
run_samtools = paste(
  'samtools view -S -b align.sam > align.bam \n
  samtools view -b -F4 align.bam > align.final.bam \n
  samtools sort align.final.bam > align.sort.bam \n
  samtools index align.sort.bam \n') #bamtofastq this file?
system(run_samtools)

# Set Path
system('git clone https://github.com/rsh249/PMDtools') #python3 compatible code?

#Deamination Mapping using PMDTools
system('gunzip *.fna.gz')
system('cat plastid.*.fna > plastid_all.fna')

run_pmdtools = paste('samtools calmd -b align.sort.bam', ref, '| samtools view -h - | python3 ../PMDtools/pmdtools.0.60.py --threshold 0 --header --deamination --range 125 --CpG > out.CpG')
system(run_pmdtools)
run_pmdtools2 = paste('samtools view -h align.sort.bam | python3 ../PMDtools/pmdtools.0.60.py --deamination --range 125 > out2.deam') #as in Moore etal. 2020
system(run_pmdtools2)

#Filtering on PMD Score of "0"
run_pmdtools_filter = paste('samtools view -h align.sort.bam | python3 ../PMDtools/pmdtools.0.60.py --threshold 0 --header | samtools view -Sb - > filter.bam') 
system(run_pmdtools_filter)

# Bam to Fastq Conversion
run_bam2fastq = paste('samtools fastq -@ 8 filter.bam > filter.fastq')
system(run_bam2fastq)

########

# Plotting (needs work)
library(ShortRead)
library(ggplot2)

raw_kreport = data.table::fread('raw_kreport.tab')
library(dplyr)
library(ggplot2)
filter_raw_kreport = raw_kreport_2 %>%
  filter(V1>=1) %>% filter(V4=='F')
ggplot(filter_raw_kreport) + geom_col(aes(x=V6, y = V2)) + theme(axis.text.x = element_text(angle=90))

filtered_kreport = data.table::fread('filtered_kreport.tab')
library(dplyr)
library(ggplot2)
filter_filtered_kreport = filtered_kreport %>%
  filter(V1>=1) %>% filter(V4=='F')
ggplot(filter_filtered_kreport) + geom_col(aes(x=V6, y = V2)) + theme(axis.text.x = element_text(angle=90))

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

#overlap reads and plot length distribution
#out2.deam
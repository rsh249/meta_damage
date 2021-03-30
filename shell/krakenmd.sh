# script to build kraken2 database

# mkdir krakendb

# kraken2-build --download-taxonomy --db krakendb
# kraken2-build --add-to-library mito_reference.fna --db krakendb
# kraken2-build --build --threads 12 --db krakendb

# kraken2 --db krakendb --threads 8 --use-names --report raw_kreport.tab --fastq-input SRR7774472.fastq > raw_kraken.out
# kraken2 --db krakendb --threads 8 --use-names --report filtered_kreport.tab --fastq-input run/filter.fastq > filtered_kraken.out

raw_kreport = data.table::fread('raw_kreport')
library(dplyr)
library(ggplot2)
filter_raw_kreport = raw_kreport %>%
   filter(V1>=0.01) %>% filter(V4=='F')
ggplot(filter_raw_kreport) +
   geom_col(aes(x=V6, y = V2) # I think this might be wrong let's try
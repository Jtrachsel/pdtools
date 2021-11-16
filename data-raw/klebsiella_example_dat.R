## code to prepare `klebsiella_example_dat` dataset goes here
library(tidyverse)
library(pdtools)
list_PDGs('Klebsiella')
system('mkdir data')
download_most_recent_complete('Klebsiella', folder_prefix = './data/')
dat <- read_tsv('./data/PDG000000012.1049.amr.metadata.tsv')

set.seed(1)
klebsiella_example_dat <-
  dat %>%
  filter(isolation_source != 'NULL') %>%
  filter(host != 'NULL') %>%
  slice_sample(n = 100)

usethis::use_data(klebsiella_example_dat, overwrite = TRUE)

#######

curl::curl_download('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt',
                    destfile = './data/assembly_summary.txt')


assembly_summary_example <- read_tsv('./data/assembly_summary.txt', skip=1) |>
  filter(`# assembly_accession` %in% klebsiella_example_dat$asm_acc)



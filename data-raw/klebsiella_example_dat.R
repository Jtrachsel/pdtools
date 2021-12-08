## code to prepare `klebsiella_example_dat` dataset goes here
library(tidyverse)
library(pdtools)
list_PDGs('Klebsiella')
system('mkdir data')
download_most_recent_complete('Klebsiella', folder_prefix = './data/')
dat <- read_tsv('./data/PDG000000012.1053.amr.metadata.tsv')
clusts <- read_tsv('./data/PDG000000012.1053.cluster_list.tsv')
dat <- dat %>% left_join(clusts)

set.seed(1)
klebsiella_example_dat <-
  dat %>%
  filter(isolation_source != 'NULL') %>%
  filter(host != 'NULL') %>%
  slice_sample(n = 200)

usethis::use_data(klebsiella_example_dat, overwrite = TRUE)

#######

curl::curl_download('https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt',
                    destfile = './data/assembly_summary.txt')


assembly_summary_example <- read_tsv('./data/assembly_summary.txt', skip=1) |>
  filter(`# assembly_accession` %in% klebsiella_example_dat$asm_acc)




#### country vector #####
library(tidyverse)
library(maps)

country_vector <-
  world.cities$country.etc %>%
  unique() %>%
  sort()

names(country_vector) <- country_vector

country_vector <-
  country_vector %>%
  sub('UK', 'UK|United Kingdom', .) %>%
  sub('USA', 'USA|United States|United States of America', .)

names(country_vector)

usethis::use_data(country_vector, overwrite = TRUE)


## code to prepare `klebsiella_example_dat` dataset goes here
library(tidyverse)
library(pdtools)

list_PDGs('Klebsiella')

usethis::use_directory('data')

download_most_recent_complete('Klebsiella', folder_prefix = './data/')
download_gbk_assembly_summary('./data/assembly_summary.txt')

files <- list.files('./data/', pattern = 'PDG', full.names = T)

dat <-
  read_tsv(files[1]) |>
  left_join(read_tsv(files[2]))

dat <-
  dat|>
  make_ftp_paths(assembly_summary_path = './data/assembly_summary.txt')



set.seed(1)

klebsiella_example_dat <-
  dat %>%
  filter(asm_acc != 'NULL') |>
  filter(isolation_source != 'NULL') %>%
  filter(host != 'NULL') %>%
  slice_sample(n = 200)

klebsiella_example_dat |> pull(ftp_path)

usethis::use_data(klebsiella_example_dat, overwrite = TRUE)

file.remove(files)
file.remove('./data/assembly_summary.txt')

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

#### gbk assem sum
# assum <-
#   read_tsv('tests/assembly_summary.txt', skip=1) |>
#   mutate(asm_acc=`# assembly_accession`) |>
#   filter(asm_acc %in% klebsiella_example_dat$asm_acc) |>
#   select(-asm_acc) |>
#   write_tsv('./tests/tmp_kleb_assembly_summary.txt', )
#
#
#
# ftps |> filter(asm_acc %in% klebsiella_example_dat$asm_acc) |>
#   write_tsv('./tests/kleb_ex_assembly')
#
# system('echo "#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file." > tests/kleb_assembly_summary.txt')
# system('cat ./tests/tmp_kleb_assembly_summary.txt >> tests/kleb_assembly_summary.txt')
# system('rm ./tests/assembly_summary.txt ./tests/tmp_kleb_assembly_summary.txt')
#

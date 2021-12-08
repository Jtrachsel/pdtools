pdtools
================

## Installation

*Requires R &gt;= 4.1* *Requires biostrings, need to install from
biocmanager*

``` r
BiocManager::install('Biostrings')
remotes::install_github('jtrachsel/pdtools')
```

## Description

A collection of functions for working with data from the [NCBI Pathogen
Detection project](https://www.ncbi.nlm.nih.gov/pathogens/)

<https://jtrachsel.github.io/pdtools/index.html>

## Examples

#### List available organisms

``` r
library(pdtools)

list_organisms()
```

#### Download the most recent metadata for an organism:

``` r
system('mkdir data')
download_most_recent_complete('Campylobacter', folder_prefix = './data/')
```

#### Join the downloaded files

``` r
# The names of these files will change based on the most recent complete data
# you will have to adjust these
meta <- readr::read_tsv('./data/PDG000000003.1540.amr.metadata.tsv') %>% 
  dplyr::left_join(readr::read_tsv('./data/PDG000000003.1540.cluster_list.tsv'))
```

#### Extract agricultural host species from several metadata fields:

``` r
# create a two column tibble containing consensus host for each isolate
host_info <- extract_consensus_ag_species(meta)

# join back to metadata
meta <- meta %>% left_join(host_info)
```

#### Extract earliest year from 3 date columns

``` r
earliest_year <- meta %>% extract_earliest_year()

meta <- meta %>% left_join(earliest_year)
```

#### Generate ftp download paths for a collection of isolates

``` r
# download most recent assembly summary
curl::curl_download('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt',
                    destfile = './data/assembly_summary.txt')

# only download isolates from swine after 2015
meta_filt <-
  meta %>%
  filter(Year > 2015 & ag_match == 'Swine') %>% 
  write_tsv('./data/swine_2015_meta.tsv')

ftp_paths <- 
  meta_filt %>% 
  make_fna_urls(.$asm_acc) %>% 
  write_lines('./data/ftp_download_paths.txt')

# download from the command line:
# cd data
# cat ftp_download_paths.txt | parallel -j 2 'wget {}'

# or within R:
# set names to the desired path 
names(ftp_paths) <-paste0('./data/',meta_filt$asm_acc, '.fna')

Map(function(u, d) download.file(u, d, mode="wb"), ftp_paths, names(ftp_paths))
```

#### Generate an input file for caclulating a pangenome with [ppanggolin](https://github.com/labgem/PPanGGOLiN)

``` r
# if you have some reference genomes that are complete (circularized) you can 
# feed their paths into the 'complete_genome_paths parameter and the function
# will correctly specify cirular contigs for ppanggolin.  

fna_files <- list.files('./data/', '.fna', full.names = T)

build_ppanggolin_file_fastas(incomplete_genome_paths = fna_files) %>% 
  write_tsv('ppanggolin_file.tsv')
```

#### Select a representative set of isolates from a pangenome

``` r
# Read in presence/absence matrix and format correctly:
pan_PA <-
  read_tsv('./pan/gene_presence_absence.Rtab')  %>% 
  column_to_rownames(var = 'Gene')  %>% 
  as.matrix() 

# this will return a a small set of genomes that contain at least the proportion of
# genes you specify

# this set of genomes will contain 99% of all the genes detected in the pangenome
get_pangenome_representatives(pan_mat = pan_PA, SEED = 2, desired_coverage = .99)
```

## TODO

-   update\_collection() function?
    -   Should take and old metadata file and a new metadata file as
        inputs.  
    -   return a vector of new genomes that were not present in the old
        list.
    -   also return genomes with newer assembly accession version  
-   extract\_consensus\_ag\_species currently assumes isolates are from
    humans if the epi\_type is clinical and no other information is
    available. This is probably wrong in some cases…
-   Reference README? <https://ftp.ncbi.nlm.nih.gov/pathogen/ReadMe.txt>

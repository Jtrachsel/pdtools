#' List available PDG accessions for an organism
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#'
#' @return returns a tibble or PDG accessions and release dates
#' @export
#'
#' @examples  #list_PDGs('Klebsiella')
#' @importFrom rlang .data
list_PDGs <- function(organism){
  #Checks the NCBI Path Det Database for the most recent version number
  # Returns a nicely formatted table
  #   Listing all available PDGs and their release dates
  # looks like the one we want is at the top...

  PDD_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/', organism)

  PDGs <- rvest::read_html(PDD_url) %>%
    rvest::html_text2() %>%
    stringr::str_split(pattern = '-PDG') %>%
    base::unlist()
  PDGs <- PDGs[-c(1, length(PDGs))]

  PDG_table <-
    tibble::tibble(raw=PDGs,  # 1st and last lines are not PDGs
            PDG=paste('PDG',sub('(.*)/(.*)','\\1',.data$raw), sep = ''),
            release_date=lubridate::ymd_hm(sub('(.*)/(.*)','\\2',.data$raw))) %>%
    dplyr::select(-.data$raw) %>%
    dplyr::arrange(dplyr::desc(.data$release_date))

  return(PDG_table)

}




#' Download Pathogen Detection metadata for a given organism
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#' @param PDG The Pathogen Detection accession number to download
#' @param folder_prefix a string to append to the download path, ie './data/'
#'
#' @return returns nothing, but probably should. Will initiate a download of the specified files
#' @export
#'
#' @examples #download_PDD_metadata(organism = 'Campylobacter',PDG = 'PDG000000003.1517')
download_PDD_metadata <- function(organism, PDG, folder_prefix=NULL){
  # browser()

  # some of these files are large, change timeout to 10 min within this function
  # the options call to change the timeout returns the original options
  original_options <- options(timeout = 600)
  on.exit(options(original_options))

  # Given an organism and a PDG accession
  # downloads the 3 metadata files from ncbi
  # seems like I dont need the "metadata" becuase the AMR table has the same info
  # meta_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism, '/',PDG,'/Metadata/',PDG,'.metadata.tsv')
  # meta_dest <- paste0(folder_prefix, PDG, '.metadata.tsv')
  options(timeout = max(500, getOption("timeout")))
  # wget https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/"$PDG"/AMR/"$PDG".amr.metadata.tsv
  amr_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/AMR/',PDG,'.amr.metadata.tsv')
  amr_dest <- paste0(folder_prefix, PDG, '.amr.metadata.tsv')

  # wget https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/"$PDG"/Clusters/"$PDG".reference_target.cluster_list.tsv
  cluster_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/Clusters/',PDG,'.reference_target.cluster_list.tsv')
  cluster_dest <- paste0(folder_prefix, PDG, '.cluster_list.tsv')

  # print('downloading metadata...')
  # curl_download(url = meta_url, destfile = meta_dest)
  print('downloading amr data...')
  utils::download.file(url = amr_url, destfile = amr_dest)
  print('downloading cluster data...')
  utils::download.file(url = cluster_url, destfile = cluster_dest)
}


#' Download the most recent complete metadata for a specified organism
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#' @param folder_prefix a string to append to the download path, ie './data/'
#'
#' @return Returns nothing, but probably should
#' @export
#'
#' @examples #download_most_recent_complete('Salmonella')
download_most_recent_complete <- function(organism, folder_prefix=NULL){

  PDG <- find_most_recent_complete(organism = organism)

  msg <- base::paste('downloading',organism, PDG[1], 'released', PDG[2])
  base::print(msg)
  download_PDD_metadata(organism = organism, PDG = PDG[1], folder_prefix = folder_prefix)

}


#' Find most recent complete PDG
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#'
#' @return returns a vector of length 2, 1=PDG accession of most recent complete, 2=release date
#' @export
#'
#' @examples #find_most_recent_complete('Salmonella')
find_most_recent_complete <- function(organism){
  # should go through the list of available PDGs and look for presence of
  # amr and cluster metadata as well
  # browser()
  url <- base::paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/', organism)
  PDG_table <- list_PDGs(organism = organism)

  index <- 1
  complete_PDG <- FALSE
  while (complete_PDG == FALSE) {
    CURRENT_PDG <- PDG_table$PDG[index]
    RELEASE_DATE <- PDG_table$release_date[index]
    complete_PDG <- check_complete_PDG(organism, PDG_table$PDG[index])
    index <- index+1
  }
  return(c(CURRENT_PDG, RELEASE_DATE))
}


#' Check if a single PDG is complete (AMR and Cluster metadata)
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#' @param PDG the ncbi pathogen detection accession to check
#'
#' @return TRUE/FALSE depending on if all URLs for download exist
#' @export
#'
#' @examples #check_complete_PDG(URL)
#'
check_complete_PDG <- function(organism, PDG){
  # browser()

  amr_url <-paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/AMR/',PDG,'.amr.metadata.tsv')
  clusters_url <-  paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/Clusters/',PDG,'.reference_target.cluster_list.tsv')

  all_urls_exist <-
    c(amr_url, clusters_url) %>%
    RCurl::url.exists() %>%
    base::all()

  return(all_urls_exist)
}


#' make specific ftp download paths for a dataframe with ftp_paths and assembly accessions
#'
#' @param type type of download path to generate, one of: 'fna', 'gbff', 'gff', 'gtf', 'faa', 'cds'
#' @param data a dataframe with the columns 'ftp_path' and 'asm_acc'
#' ftp_path should be a column produced by the function make_ftp_paths()
#'
#' @return returns the original dataframe with an added column, named "{type}_download"
#' @export
#'
#' @examples #make_download_urls(asm_acc = klebsiella_example_dat$asm_acc,
#'                              ftp_paths = klebsiella_example_dat$ftp_path, type='fna')
#'


make_download_urls <- function(data, type){

  suffixes=base::c(fna='_genomic.fna.gz',
                   gbff='_genomic.gbff.gz',
                   gff='_genomic.gff.gz',
                   gtf='_genomic.gtf.gz ',
                   faa='_protein.faa.gz',
                   cds='_cds_from_genomic.fna.gz')

  if (!(type %in% base::names(suffixes))){
    base::errorCondition(base::paste0('"type" must be one of ','"', base::paste(base::names(suffixes), collapse = ' '),'"'))
  }



  result <-
    data %>%
    mutate("{type}_download":=
             base::paste0(ftp_paths,
                          '/',
                          base::sub('https://ftp.ncbi.nlm.nih.gov/genomes/all/.*/.*/.*/.*/(.*)', '\\1',
                                    ftp_paths),
                          suffixes[type]))
  return(result)

}


#' generate ftp site download urls for all SNP trees containing the provided isolates
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#' @param data a metadata table, must contain the column 'PDS_acc' from merging in the SNP cluster data
#' @param PDG The PDG version the metadata is from.
#'
#' @return returns a vector of ftp download urls for each tar.gz file containing the SNP tree info
#' @export
#'
#' @examples make_SNPtree_urls(organism = 'Klebsiella',
#'  data = klebsiella_example_dat, PDG = 'PDG000000012.1053')
make_SNPtree_urls <- function(organism, data, PDG){
  # One SNP tree for each PDG represented in the data
  # Organism <- 'Klebsiella'
  # PDG <- 'PDG000000012.1053'

  num_no_clust <- base::sum(base::is.na(data$PDS_acc))

  PDSs <- data %>% dplyr::filter(!is.na(.data$PDS_acc)) %>% dplyr::pull(.data$PDS_acc) %>% base::unique()
  urls <- base::paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/', PDG, '/SNP_trees/', PDSs, '.tar.gz')
  base::message(base::paste(num_no_clust, 'Isolates in the collection are not represented in SNP trees'))
  return(urls)
}

#' List organisms available in the Pathogens database
#'
#' @return a tibble of available organisms
#' @export
#'
#' @examples  list_organisms()
list_organisms <- function(){
  url <- 'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/'

  organisms <-
    rvest::read_html(url) %>%
    rvest::html_text2() %>%
    stringr::str_split(pattern = ' -') %>%
    base::unlist()
  organisms <- organisms[-c(1, length(organisms))]

  organism_table <-
    tibble::tibble(raw=organisms,
                   organism=base::sub('(.*)/(.*)','\\1',.data$raw),
                   release_date=lubridate::ymd_hm(sub('(.*)/(.*)','\\2',.data$raw))) %>%
    dplyr::select(-.data$raw)

  return(organism_table)


}



#' generate ftp site paths for a selection of assembly accessions
#'
#' @param asm_accessions vector of assembly accessions
#' @param assembly_summary_path path to genbank assembly_summary.txt, see download_gbk_assembly_summary()
#'
#'
#' @return a two column tibble 1= asm_acc ; 2= ftp_path
#' @export
#'
#' @examples #make_ftp_paths(klebsiella_example_data$asm_acc, './test/kleb_assembly_summary.txt')
make_ftp_paths <- function(data, assembly_summary_path){

  # browser()

  ftp_asm_map <-
    readr::read_tsv(assembly_summary_path, skip=1) %>%
    dplyr::transmute(asm_acc=.data$`# assembly_accession`,
                     .data$ftp_path)
  data %>% left_join(ftp_asm_map)

  return(result)

}

#' Convenience function to download the assembly_summary.txt file from genbank
#'
#' @param destfile passed to download.file()'s destfile, path to store the downloaded file
#'
#' @return returns nothing but probably should...
#' @export
#'
#' @examples #not run download_gbk_assembly_summary(destfile='assembly_summary.txt')
download_gbk_assembly_summary <- function(destfile){
  original_options <- base::options(timeout = 6000)
  base::on.exit(base::options(original_options))

  utils::download.file('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt',
                      destfile = destfile)

}



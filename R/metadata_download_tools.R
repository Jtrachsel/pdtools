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
  original_options <- options(timeout = 6000)
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


  if(!file.exists(amr_dest)){
    print('downloading amr data...')
    utils::download.file(url = amr_url, destfile = amr_dest)
  } else {
    print(paste('file', amr_dest, 'already exists...', 'skipping'))
  }

  if(!file.exists(cluster_dest)){
    print('downloading cluster data...')
    utils::download.file(url = cluster_url, destfile = cluster_dest)
  }else {
    print(paste('file', cluster_dest, 'already exists...', 'skipping'))
  }


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
#' @noRd
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
#' @noRd
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



######## NEW! ###########

#' get the version of the PDD metadata
#'
#' @param data_dir directory that contains the downloaded PDD metadata
#'
#' @return character vector length 1, PDG accession
#'
#'
#' @examples # get_PDG_version('./data/')
get_PDG_version <- function(data_dir){
  base::sub('(PDG.*).amr.metadata.tsv','\\1',base::list.files(data_dir, '(PDG.*).amr.metadata.tsv')[1])
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



#' Download SNP trees from a dataframe containing SNP_tree_urls and dests
#'
#' @param data dataframe with SNP tree URLs and dests
#'
#' @return the input dataframe with an added column indicating the status of the download
#' @export
#'
#' @examples # soon
download_SNP_trees <- function(data){

  original_options <- base::options(timeout = 10000)
  base::on.exit(base::options(original_options))

  ### check if SNPS exist here
  data %>%
    dplyr::mutate(SNP_tree_dl=purrr::map2(.x = .data$SNP_tree_url, .y=.data$SNP_tree_dest, .f = ~utils::download.file(.x, .y)))
}
#
# make_SNPtree_urls(organism = 'Salmonella',
#                   data = infmeta,
#                   PDG = get_PDG_version('metadata'))

#' Make SNP tree download destination paths
#'
#' @param data A dataframe containing the download paths for desired SNP trees
#' @param data_dir a path to a directory to store the downloaded trees in
#'
#' @return the input dataframe with an added column containing the SNP_tree_dest column
#' @export
#'
#' @examples # soon
make_SNP_tree_dest <- function(data, data_dir){
  data %>%
    dplyr::mutate(SNP_tree_dest=
             base::paste0(data_dir, '/', base::sub('.*SNP_trees/(PDS[0-9]+.[0-9]+.tar.gz)','\\1',.data$SNP_tree_url)))
}



# ADD function to download most up to date AMR reference gene catalogue
# need to look for most recent version 1st
#  https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.10/2021-12-21.1/ReferenceGeneCatalog.txt



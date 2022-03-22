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



######## NEW! ###########

get_PDG_version <- function(data_dir){
  sub('(PDG.*).amr.metadata.tsv','\\1',list.files(data_dir, '(PDG.*).amr.metadata.tsv')[1])
}


download_SNP_trees <- function(data){

  original_options <- base::options(timeout = 10000)
  base::on.exit(base::options(original_options))

  ### check if SNPS exist here
  data %>%
    mutate(SNP_tree_dl=map2(.x = SNP_tree_url, .y=SNP_tree_dest, .f = ~download.file(.x, .y)))
}
#
# make_SNPtree_urls(organism = 'Salmonella',
#                   data = infmeta,
#                   PDG = get_PDG_version('metadata'))

make_SNP_tree_dest <- function(data, data_dir){
  data %>%
    mutate(SNP_tree_dest=
             paste0(data_dir, '/', sub('.*SNP_trees/(PDS[0-9]+.[0-9]+.tar.gz)','\\1',SNP_tree_url)))
}




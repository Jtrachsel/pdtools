#' List available PDG accessions for an organism
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#'
#' @return returns a tibble or PDG accessions and release dates
#' @export
#'
#' @examples  list_PDGs('Klebsiella')
list_PDGs <- function(organism){
  #Checks the NCBI Path Det Database for the most recent version number
  # Returns a nicely formatted table
  #   Listing all available PDGs and their release dates
  # looks like the one we want is at the top...

  PDD_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/', organism)

  PDGs <- read_html(PDD_url) %>%
    rvest::html_text2() %>%
    str_split(pattern = '-PDG') %>%
    unlist()
  PDGs <- PDGs[-c(1, length(PDGs))]

  PDG_table <-
    tibble(raw=PDGs,  # 1st and last lines are not PDGs
           PDG=paste('PDG',sub('(.*)/(.*)','\\1',raw), sep = ''),
           release_date=ymd_hm(sub('(.*)/(.*)','\\2',raw))) %>%
    select(-raw) %>%
    arrange(desc(release_date))

  return(PDG_table)

}




#' Download Pathogen Detection metadata for a given organism
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#' @param PDG The Pathogen Detection accession number to download
#' @param folder_prefix a string to append to the download path, ie './data/'
#'
#' @return
#' @export
#'
#' @examples download_PDD_metadata(organism = 'Campylobacter',PDG = 'PDG000000003.1517')
download_PDD_metadata <- function(organism, PDG, folder_prefix=NULL){

  # Given an organism and a PDG accession
  # downloads the 3 metadata files from ncbi
  # must have a ./data/ folder in the working directory

  meta_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism, '/',PDG,'/Metadata/',PDG,'.metadata.tsv')
  meta_dest <- paste0(folder_prefix, PDG, '.metadata.tsv')

  # wget https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/"$PDG"/AMR/"$PDG".amr.metadata.tsv
  amr_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/AMR/',PDG,'.amr.metadata.tsv')
  amr_dest <- paste0(folder_prefix, PDG, '.amr.metadata.tsv')

  # wget https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/"$PDG"/Clusters/"$PDG".reference_target.cluster_list.tsv
  cluster_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/Clusters/',PDG,'.reference_target.cluster_list.tsv')
  cluster_dest <- paste0(folder_prefix, PDG, '.cluster_list.tsv')

  print('downloading metadata...')
  curl_download(url = meta_url, destfile = meta_dest)
  print('downloading amr data...')
  curl_download(url = amr_url, destfile = amr_dest)
  print('downloading cluster data...')
  curl_download(url = cluster_url, destfile = cluster_dest)
}



#' Download the most recent complete metadata for a specified organism
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#' @param folder_prefix a string to append to the download path, ie './data/'
#'
#' @return
#' @export
#'
#' @examples download_most_recent_complete('Salmonella')
download_most_recent_complete <- function(organism, folder_prefix=NULL){

  PDG <- find_most_recent_complete(organism = organism)

  msg <- paste('downloading',organism, PDG, 'released', second_most_recent$release_date)
  print(msg)
  download_PDD_metadata(organism = organism, PDG = PDG, folder_prefix = folder_prefix)

}


#' Find most recent complete PDG
#'
#' @param organism
#'
#' @return
#' @export
#'
#' @examples find_most_recent_complete('Salmonella')
find_most_recent_complete <- function(organism){
  # should go through the list of available PDGs and look for presence of
  # amr and cluster metadata as well
  # browser()
  url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/', organism)
  PDG_table <- list_PDGs(organism = organism) %>%
    mutate(URL=paste0(url, '/', PDG))

  index <- 1
  FOUND_IT <- FALSE
  while (FOUND_IT == FALSE) {
    FOUND_IT <- check_complete_PDG(PDG_table$URL[index])
    CURRENT_PDG <- PDG_table$PDG[index]
    index <- index+1
  }
  return(CURRENT_PDG)
}


#' Check if a single PDG is complete (AMR and Cluster metadata)
#'
#' @param url
#'
#' @return
#' @export
#'
#' @examples check_complete_PDG(URL)
check_complete_PDG <- function(url){
  contents <- read_html(url) %>%
    rvest::html_text2()
  AMR <- grepl('AMR', contents)
  META <- grepl('Metadata', contents)
  CLUSTERS <- grepl('Clusters', contents)
  return(all(c(AMR, META, CLUSTERS)))
}


# find_most_recent_complete('Escherichia_coli_Shigella')



#
#
# # test on Sal and Campy
# Sal_PDGs <- list_PDGs('Salmonella')
#
#
# PDG <- Sal_PDGs %>% slice_head() %>% pull(PDG)
#
#
# download_PDD_metadata('Salmonella', PDG = PDG)
#
#
# Campy_PDGs <- list_PDGs('Campylobacter')
#
#
#
# Campy_PDG <- Campy_PDGs %>% slice_head() %>% pull(PDG)
#
# download_PDD_metadata('Campylobacter', PDG = Campy_PDG)
#
#

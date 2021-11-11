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

  PDGs <- rvest::read_html(PDD_url) |>
    rvest::html_text2() |>
    str_split(pattern = '-PDG') |>
    unlist()
  PDGs <- PDGs[-c(1, length(PDGs))]

  PDG_table <-
    tibble::tibble(raw=PDGs,  # 1st and last lines are not PDGs
            PDG=paste('PDG',sub('(.*)/(.*)','\\1',raw), sep = ''),
            release_date=lubridate::ymd_hm(sub('(.*)/(.*)','\\2',raw))) |>
    dplyr::select(-raw) |>
    dplyr::arrange(desc(release_date))

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
  # Given an organism and a PDG accession
  # downloads the 3 metadata files from ncbi
  # seems like I dont need the "metadata" becuase the AMR table has the same info
  # meta_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism, '/',PDG,'/Metadata/',PDG,'.metadata.tsv')
  # meta_dest <- paste0(folder_prefix, PDG, '.metadata.tsv')

  # wget https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/"$PDG"/AMR/"$PDG".amr.metadata.tsv
  amr_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/AMR/',PDG,'.amr.metadata.tsv')
  amr_dest <- paste0(folder_prefix, PDG, '.amr.metadata.tsv')

  # wget https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/"$PDG"/Clusters/"$PDG".reference_target.cluster_list.tsv
  cluster_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/Clusters/',PDG,'.reference_target.cluster_list.tsv')
  cluster_dest <- paste0(folder_prefix, PDG, '.cluster_list.tsv')

  # print('downloading metadata...')
  # curl_download(url = meta_url, destfile = meta_dest)
  print('downloading amr data...')
  curl::curl_download(url = amr_url, destfile = amr_dest)
  print('downloading cluster data...')
  curl::curl_download(url = cluster_url, destfile = cluster_dest)
}


#' Download the most recent complete metadata for a specified organism
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#' @param folder_prefix a string to append to the download path, ie './data/'
#'
#' @return
#' @export
#'
#' @examples #download_most_recent_complete('Salmonella')
download_most_recent_complete <- function(organism, folder_prefix=NULL){

  PDG <- find_most_recent_complete(organism = organism)

  msg <- paste('downloading',organism, PDG[1], 'released', PDG[2])
  print(msg)
  download_PDD_metadata(organism = organism, PDG = PDG[1], folder_prefix = folder_prefix)

}


#' Find most recent complete PDG
#'
#' @param organism a string ie 'Salmonella' or 'Campylobacter' etc
#'
#' @return
#' @export
#'
#' @examples #find_most_recent_complete('Salmonella')
find_most_recent_complete <- function(organism){
  # should go through the list of available PDGs and look for presence of
  # amr and cluster metadata as well
  # browser()
  url <- paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/', organism)
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
#' @param PDG
#'
#' @return TRUE/FALSE depending on if all URLs for download exist
#' @export
#'
#' @examples #check_complete_PDG(URL)
check_complete_PDG <- function(organism, PDG){
  # browser()

  amr_url <-paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/AMR/',PDG,'.amr.metadata.tsv')
  clusters_url <-  paste0('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/',organism,'/',PDG,'/Clusters/',PDG,'.reference_target.cluster_list.tsv')

  all_urls_exist <-
    c(amr_url, clusters_url) |>
    RCurl::url.exists() |>
    all()

  return(all_urls_exist)
}





#' Generate ftp urls for fna files from a vector of genbank assembly accessions
#'
#' @param asm_accessions character vector of genbank accessions
#' @param assembly_summary assembly summary downloaded from ncbi, must have asm_acc column
#'
#'
#'
#' @return returns a vector of ftp urls
#' @export
#'
#' @examples #soon
make_fna_urls <- function(asm_accessions, assembly_summary){
  # assembly summary needs to have the accessions changed from
  # `# assembly_accession` to asm_acc
  ass_sum <- ass_sum |>
    curl::filter(asm_acc %in% asm_accessions)
  paste0(ass_sum$ftp_path,
         '/',
         sub('https://ftp.ncbi.nlm.nih.gov/genomes/all/.*/.*/.*/.*/(.*)', '\\1',
             ass_sum$ftp_path),
         '_genomic.fna.gz')

}

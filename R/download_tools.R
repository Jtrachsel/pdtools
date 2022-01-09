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
  base::download.file(url = amr_url, destfile = amr_dest)
  print('downloading cluster data...')
  base::download.file(url = cluster_url, destfile = cluster_dest)
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


#' Generate a two column tibble mapping ftp paths to assembly accessions
#'
#' @param filename path to the ncbi genbank assembly_summary.txt
#'     curl -O https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
#'
#' @return returns a two column tibble of asm_acc and ftp
#'
#' @examples #ftp_paths <- ftp_paths_from_assem_sum('./data/assembly_summary.txt')
#' @importFrom rlang .data
#' @noRd
ftp_paths_from_assem_sum <- function(filename){
  # browser()
  dat <- readr::read_tsv(filename, skip=1) %>%
    dplyr::transmute(asm_acc=.data$`# assembly_accession`,
              .data$ftp_path)
  return(dat)
}



#' generate fasta ftp download paths for a vector of assembly accessions
#'
#' @param asm_accessions vector of assembly accessions to generate ftp paths for
#' @param assembly_summary_path path to an assembly_summary.txt file downloaded from ncbi
#'
#' @return returns a vector of ftp paths for use with curl wget etc.
#' @export
#'
#' @examples # tst <- make_fna_urls(klebsiella_example_dat$asm_acc, './data/assembly_summary.txt')
make_fna_urls <- function(asm_accessions, assembly_summary_path){
  # assembly summary needs to have the accessions changed from
  # `# assembly_accession` to asm_acc
  selected_ftp_paths <-
    ftp_paths_from_assem_sum(assembly_summary_path) %>%
    dplyr::filter(.data$asm_acc %in% asm_accessions)

  base::paste0(selected_ftp_paths$ftp_path,
               '/',
               base::sub('https://ftp.ncbi.nlm.nih.gov/genomes/all/.*/.*/.*/.*/(.*)', '\\1',
                         selected_ftp_paths$ftp_path),
               '_genomic.fna.gz')

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




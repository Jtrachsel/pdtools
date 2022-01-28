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


#' generate ftp site paths for a selection of assembly accessions
#'
#' @param assembly_summary_path path to genbank assembly_summary.txt, see download_gbk_assembly_summary()
#'
#' @param data a dataframe containing an asm_acc
#'
#' @return a two column tibble 1= asm_acc ; 2= ftp_path
#' @export
#'
#' @examples #make_ftp_paths(klebsiella_example_data './test/assembly_summary.txt')
make_ftp_paths <- function(data, assembly_summary_path){

  # browser()
  # should check for NAs or weirdly formatted asm_acc

  # check_asm_acc

  ftp_asm_map <-
    readr::read_tsv(assembly_summary_path, skip=1) %>%
    dplyr::transmute(asm_acc=.data$`# assembly_accession`,
                     .data$ftp_path)
  dplyr::filter(grepl('https://ftp.ncbi.nlm.nih.gov',.data$ftp_path))

  result <- data %>% dplyr::left_join(ftp_asm_map)

  return(result)

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
#' @examples make_download_urls(klebsiella_example_dat, type='fna')
#' @importFrom rlang :=
make_download_urls <- function(data, type){
  suffixes <- supported_download_types(type)
  result <-
    data %>%
    dplyr::mutate("{type}_download":=
             base::paste0(.data$ftp_path,
                          '/',
                          base::sub('https://ftp.ncbi.nlm.nih.gov/genomes/all/.*/.*/.*/.*/(.*)', '\\1',
                                    .data$ftp_path),
                          suffixes[type]))
  return(result)

}

#' Make download destination paths
#'
#' @param data A dataframe containing an asm_acc column
#' @param type the type of files you want to download, one of: 'fna', 'gbff', 'gff', 'gtf', 'faa', 'cds'
#' @param dest_dir path to the directory you want to use, must exist, should include a trailing '/'
#'
#' @return returns a dataframe with an added "{type}_dest" column containing the paths to pass to download.file
#' @export
#' @importFrom rlang :=
#'
#' @examples # download_data %>% make_dest_faths(type='fna', dest_dir='./data/')
make_dest_paths <- function(data, type, dest_dir){

  base::file.exists(dest_dir)
  supported_download_types(type)

  data %>%
    dplyr::mutate("{type}_dest":=paste0(dest_dir, .data$asm_acc, '.', type, '.gz'))
}

#' Download specified files from NCBI ftp site
#'
#' @param data a dataframe with columns created by make_download_urls() and make_dest_paths()
#' @param type the type of files you want to download, one of: 'fna', 'gbff', 'gff', 'gtf', 'faa', 'cds'
#'
#' @return the results of attempting to download the specified files,
#'  A dataframe with two added columns:
#'    1. return of download.file, should be 0
#'    2. error column from purrr::safely(), should contain any error messages.
#' @export
#'
#' @importFrom rlang :=
#' @examples # download_data %>% download_genomes('fna')
download_genomes <-
  function(data, type){
    supported_download_types(type)
    url_var <- base::paste0(type, '_download')
    dest_var <- base::paste0(type, '_dest')
    err_var <- stats::setNames(base::list(base::as.character), glue::glue("{type}_dl_error"))

    safe_download <- purrr::safely(utils::download.file)

    data %>%
      dplyr::select(.data$asm_acc, dplyr::starts_with(type)) %>%
      dplyr::mutate("{type}_dl":=purrr::map2(.x=!!rlang::sym(url_var), .y=!!rlang::sym(dest_var), .f = ~safe_download(.x, .y))) %>%
      tidyr::unnest_wider(glue::glue('{type}_dl'),
                          names_sep = '_',
                          simplify = TRUE,
                          transform = err_var)

  }








#' helper function to check if user supplied type is in the supported types
#'
#' @param type a user input string to check
#'
#' @return returns a named vector of acceptable files and they appropriate suffixes
#'
#' @examples # supported_download_types('fna')
supported_download_types <-
  function(type){
    suffixes=base::c(fna='_genomic.fna.gz',
                     gbff='_genomic.gbff.gz',
                     gff='_genomic.gff.gz',
                     gtf='_genomic.gtf.gz ',
                     faa='_protein.faa.gz',
                     cds='_cds_from_genomic.fna.gz')

    if (!(type %in% base::names(suffixes))){
      base::errorCondition(base::paste0('"type" must be one of ','"', base::paste(base::names(suffixes), collapse = ' '),'"'))
    }
    return(suffixes)
  }



#' Convenience function to download the assembly_summary.txt file from genbank
#'
#' @param destfile passed to download.file()'s destfile, path to store the downloaded file
#' @param organism limit the gbk assembly summary to only this organism
#'
#' @return result of download.file
#' @export
#'
#' @examples #not run download_gbk_assembly_summary(destfile='assembly_summary.txt')
download_gbk_assembly_summary <- function(destfile, organism=NULL){
  original_options <- base::options(timeout = 6000)
  base::on.exit(base::options(original_options))
  if(is.null(organism)){
    utils::download.file('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt',
                         destfile = destfile)
  } else {
    download_url <- base::paste0('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/', organism, '/assembly_summary.txt')
    utils::download.file(download_url, destfile = destfile)
  }



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
                     .data$ftp_path) %>%
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
  # browser()
  if(!all(grepl('https://ftp.ncbi.nlm.nih.gov/genomes/all/', data$ftp_path))){
    return(errorCondition('some of the ftp_paths are invalid, must match "https://ftp.ncbi.nlm.nih.gov/genomes/all/" '))
  }

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
  stopifnot({
    base::file.exists(dest_dir)

  })
  supported_download_types(type)
  data %>%
    dplyr::mutate("{type}_dest":=paste0(dest_dir, .data$asm_acc, '.', type, '.gz'))
}



#' Download a set of specified genome files
#'
#' @param data dataframe that both make_download_urls and make_dest_paths have run on
#' @param type one of the accepted types of files to download
#' @param PARALLEL should downloads be run in parallel, requires you to set your
#' future::plan()
#'
#' @return the original dataframe with columns added to represent the result of the download attempt
#' @export
#'
#' @examples #
download_genomes <-
  function(data, type, PARALLEL=FALSE){
    # browser()
    supported_download_types(type)

    url_var <- base::paste0(type, '_download')
    dest_var <- base::paste0(type, '_dest')
    exists_var <-base::paste0(type, "_exists")

    # this may be unecessary
    safe_download <- purrr::possibly(utils::download.file, otherwise = 1)

    # check for existing files #
    # with and without extension in case they've been unzipped
    # print('checking for existing files')
    data <- check_if_files_exist(data, type)
    need_to_download <- data %>% dplyr::filter(!(!!rlang::sym(exists_var)))
    already_existing <- data %>%
      dplyr::filter((!!rlang::sym(exists_var))) %>%
      dplyr::mutate("{type}_dl":='exists')

    if (base::nrow(need_to_download) > 0){

      base::print('downloading genomes, please be patient')

      if (PARALLEL){
        print("using parallel downloads, if you haven't set your future::plan()
              these will still be sequential")
        res <-
          need_to_download %>%
          dplyr::mutate("{type}_dl":=furrr::future_map2_chr(.x=!!rlang::sym(url_var), .y=!!rlang::sym(dest_var), .f = ~safe_download(url=.x, destfile=.y, quiet=TRUE)))

        res <- dplyr::bind_rows(already_existing, res) %>% check_if_files_exist(type)

        return(res)

      } else {
        res <- need_to_download %>%
          dplyr::mutate("{type}_dl":=purrr::map2_chr(.x=!!rlang::sym(url_var), .y=!!rlang::sym(dest_var), .f = ~safe_download(url=.x, destfile=.y, quiet=TRUE))) #%>%
        res <- dplyr::bind_rows(already_existing, res)

        return(res)
      }


    } else {

      return(already_existing)

    }

  }






#' helper function to check if user supplied type is in the supported types
#'
#' @param type a user input string to check
#'
#' @return returns a named vector of acceptable files and they appropriate suffixes
#' @noRd
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
      return(base::errorCondition(base::paste0('"type" must be one of ','"', base::paste(base::names(suffixes), collapse = ' '),'"')))
    }
    return(suffixes)
  }

download_reference_genomes <- function(genome_names,type, data_dir, avail=FALSE){

  references <- base::c(
    LT2='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic',
    Enteriditis='ENTERIDITIS',
    USDA15WA1='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/874/805/GCA_006874805.1_ASM687480v1/GCA_006874805.1_ASM687480v1_genomic',
    FSIS1502916='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/931/575/GCF_001931575.1_ASM193157v1/GCF_001931575.1_ASM193157v1_genomic')
  if (avail){
    return(references)
  }
  # browser()
  dl_tib <-
   tibble::tibble(NAME=genome_names,
           download_path=base::paste0(references[.data$NAME],'.', type, '.gz'),
           dest_path=base::paste0(data_dir,'/', .data$NAME,'.', type,'.gz')) %>%
    dplyr::mutate(RESULT=purrr::map2_int(.x=.data$download_path, .y=.data$dest_path, .f=~utils::download.file(.x,.y)))
  return(dl_tib)


}


check_if_files_exist <- function(data, type){
  dest_var <- base::paste0(type, '_dest')
  # err_var <- stats::setNames(base::list(base::as.character), glue::glue("{type}_dl_error"))
  exists_var <- glue::glue("{type}_exists")

  # check for existing files #
  # with and without extension in case they've been unzipped
  print('checking for existing files')
  exist_dat <-
    data %>%
    dplyr::select(.data$asm_acc, dest_var) %>%
    dplyr::mutate(gunzipped=base::sub('.gz','', !!rlang::sym(dest_var))) %>%
    dplyr::mutate(EXISTS=purrr::map_lgl(.x = !!rlang::sym(dest_var), .f = base::file.exists),
                  EXISTS2=purrr::map_lgl(.x = .data$gunzipped, .f = base::file.exists),
                  "{type}_exists":=.data$EXISTS | .data$EXISTS2) %>%
    dplyr::select(.data$asm_acc, dplyr::all_of(dest_var), dplyr::all_of(exists_var))

  data <- data %>% dplyr::left_join(exist_dat) %>% unique()

  if (base::all(exist_dat[[exists_var]])){
    base::print(glue::glue('all of the {type} files exist'))
  }

  if (base::any(exist_dat[[exists_var]])){
    base::print(glue::glue('some of the {type} files already exist, skipping these'))
  }

  return(data)
}


# after downloading files, make sure no asm_acc is repeated
# if repeats, keep newest accession
resolve_updated_assemblies <- function(data_dir){
  list.files(data_dir)
}



#' Return a string of names of the matching regex patterns
#'
#' @param pattern_vec a named vector of regex patterns to use
#' @param search_string a sting to search for matches
#'
#' @return names of matching patterns concatenated with '_'
#'
#'
#' @examples #pat_vec <- c(Swine='hog|swine|sow'); pdtools:::return_ag_match(pattern_vec=pat_vec, 'Hog')
#' @noRd
matches_from_vector_of_patterns <-

  function(pattern_vec, search_string){
    # was return_ag_match()
    # browser()

  match_vec <- purrr::map_lgl(.x =pattern_vec, .f=~base::grepl(.x, search_string, ignore.case = TRUE))

  res <- names(match_vec)[match_vec] %>%
    base::paste(collapse = '_')
  return(res)
}


#' Extract a consensus ag host species from metadata
#' This is pretty crude and should not be relied on for 100% accuracy
#'
#'
#' @param dat an ncbi pathogen detection metadata table
#' @param parallel boolean, should furrr be used to parallelize? need to set your future::plan()
#'
#' @return returns a tibble of 2 columns, 1st = target_acc, 2nd = ag_match
#' @export
#'
#' @examples extract_consensus_ag_species(klebsiella_example_dat)
#' @importFrom rlang .data
extract_consensus_ag_species <- function(dat, parallel=FALSE){
  # browser()
  pattern_vec <-
    base::c(Swine="swine|pork|porcine|sow|sus|hog|pig|scrofa|\\bham\\b",
      Bovine="bovine|beef|veal|cow|cattle|bos|steer|taurus|calf|bull",
      Chicken="chicken|chick|gallus|broiler|egg",
      Turkey="turkey|meleagris|gallopavo",
      Human="human|homo|sapiens|infant",
      Horse="equine|equus|horse|caballus",
      Dog="canine|dog|canis|canine",
      Cat= "\\bcat\\b|felis|catus|feline",
      Goat="caprine|goat",
      Sheep="\\bovine|sheep|lamb",
      Duck="duck",
      Goose='goose')

  # check for needed columns #
  stopifnot({

    check_for_columns(dat, c('isolation_source', 'host', 'ontological_term', 'epi_type'))

  })

  #

  if (parallel){
    first_pass <-
      dat %>%
      dplyr::transmute(target_acc=.data$target_acc,
                       search_vals=base::paste(.data$isolation_source, .data$host, .data$ontological_term,.data$epi_type)) %>%
      dplyr::mutate(ag_match=furrr::future_map_chr(.x = .data$search_vals, ~matches_from_vector_of_patterns(pattern_vec, search_string = .x)))


  } else {

    first_pass <-
      dat %>%
      dplyr::transmute(target_acc=.data$target_acc,
                       search_vals=base::paste(.data$isolation_source, .data$host, .data$ontological_term,.data$epi_type)) %>%
      dplyr::mutate(ag_match=purrr::map_chr(.x = .data$search_vals, ~matches_from_vector_of_patterns(pattern_vec, search_string = .x)))


  }

    finished <- first_pass %>%
      dplyr::filter(.data$ag_match != '')

    if (base::nrow(finished) == base::nrow(dat)) {
      result <- finished %>% dplyr::select('target_acc', 'ag_match')
      return(result)
    } else {
      # This assumes that any clinical isolate without other info is a human isolate
      second_pass <- first_pass %>%
        dplyr::filter(.data$ag_match == '') %>%
        dplyr::mutate(ag_match=base::ifelse(base::grepl('clinical', .data$search_vals), 'Human', 'Other'))


      result <- dplyr::bind_rows(finished, second_pass) %>%
        dplyr::select('target_acc', 'ag_match')

      return(result)

    }




}

#' check that all needed columns exist
#'
#' @param data dataframe to check
#' @param column_names which column names should exist?
#'
#' @return TRUE/FALSE
#' @noRd
#'
#' @examples # soon
check_for_columns <- function(data, column_names){
    all(column_names %in% colnames(data))

}
#' return a Year column containing the earliest year from the
#' available 'date' fields.
#'
#' @param PDD_metadata_table an ncbi pathogen detection metadata table
#'
#' @return returns the a tibble with a `Year` column,
#' @export
#'
#' @examples extract_earliest_year(klebsiella_example_dat)
#' @importFrom rlang .data
extract_earliest_year <- function(PDD_metadata_table){
  # was return_earliest_year
    result <-
    PDD_metadata_table %>%
    dplyr::select('target_acc', dplyr::ends_with('date')) %>%
    dplyr::mutate(dplyr::across(.cols = dplyr::ends_with('date'), .fns = base::as.character)) %>%
    tidyr::pivot_longer(cols = dplyr::ends_with('date'), names_to = 'type', values_to = 'date') %>%
    dplyr::mutate(year = base::as.numeric(base::sub('([0-9][0-9][0-9][0-9]).*','\\1',.data$date))) %>%
    dplyr::group_by(.data$target_acc) %>%
    dplyr::summarise(Year=.data$year[base::which.min(.data$year)])

  return(result)

}




#' extract collecting agency from several metadata fields
#'
#' @param meta an ncbi pathogen detection metadata table
#' @param parallel boolean, should furrr be used to parallelize? need to set your future::plan()
#'
#' @return a tibble with 2 columns, 1 = target_acc, 2 = collection_agency
#' @export
#'
#' @examples klebsiella_example_dat %>% extract_collection_agency()
extract_collection_agency <-
  function(meta, parallel=FALSE){
  pattern_vec <-
    base::c(CDC="CDC|Center for Disease Control",
            FDA="FDA|Food and Drug Administration",
            FSIS="FSIS|Food Saftey Inspection Service",
            USDA="USDA|United States? Department of Agriculture")



  if (parallel){
    first_pass <-
      meta %>%
      dplyr::transmute(target_acc=.data$target_acc,
                       search_vals=base::paste(.data$collected_by, .data$bioproject_center, sep = '_')) %>%
      dplyr::mutate(collection_agency=furrr::future_map_chr(.x = .data$search_vals, ~matches_from_vector_of_patterns(pattern_vec, search_string = .x)))

  } else {
    first_pass <-
      meta %>%
      dplyr::transmute(target_acc=.data$target_acc,
                       search_vals=base::paste(.data$collected_by, .data$bioproject_center, sep = '_')) %>%
      dplyr::mutate(collection_agency=purrr::map_chr(.x = .data$search_vals, ~matches_from_vector_of_patterns(pattern_vec, search_string = .x)))

  }


  finished <- first_pass %>%
    dplyr::filter(.data$collection_agency != '') %>%
    dplyr::select('target_acc', 'collection_agency')
  return(finished)
}

#' Extract a standardized country name from geo_loc_tag column
#'
#' @param meta A PDD metadata table
#' @param parallel boolean, should furrr be used to parallelize? need to set your future::plan()
#'
#' @return a two column tibble with target_acc and country as columns
#' @export
#'
#' @examples klebsiella_example_dat %>% extract_country()
extract_country <- function(meta, parallel=FALSE){
  # return an acceptable country column from a metadata table

  pattern_vec <- pdtools::country_vector

  if (parallel){
    first_pass <-
      meta %>%
      dplyr::transmute(target_acc=.data$target_acc,
                       search_vals=base::tolower(.data$geo_loc_name)) %>%
      dplyr::mutate(country=furrr::future_map_chr(.x = .data$search_vals, ~matches_from_vector_of_patterns(pattern_vec, search_string = .x)))

  } else {
    first_pass <-
      meta %>%
      dplyr::transmute(target_acc=.data$target_acc,
                       search_vals=base::tolower(.data$geo_loc_name)) %>%
      dplyr::mutate(country=purrr::map_chr(.x = .data$search_vals, ~matches_from_vector_of_patterns(pattern_vec, search_string = .x)))

  }

  finished <- first_pass %>%
    dplyr::filter(.data$country != '') %>%
    dplyr::select('target_acc', 'country')
  return(finished)

}

#### WORK ON THIS

#' Extract the state name from the geo_loc_name column, this only makes sence for USA isolates
#'
#' @param data a dataframe containing a `geo_loc_name` column
#'
#' @return a tibble with 2 columns, `target_acc` and `State`
#' @export
#'
#' @examples klebsiella_example_dat %>% extract_state()
extract_state <- function(data){

  state_vector <- base::paste('\\b',datasets::state.name,'\\b', '|', '\\b', datasets::state.abb,'\\b', sep = '')
  base::names(state_vector) <- datasets::state.abb
  pattern_vec <- state_vector

  first_pass <-
    data %>%
    dplyr::transmute(target_acc=.data$target_acc,
                     search_vals=base::tolower(.data$geo_loc_name)) %>%
    dplyr::mutate(State=furrr::future_map_chr(.x = .data$search_vals, ~matches_from_vector_of_patterns(pattern_vec, search_string = .x))) %>%
    dplyr::select('target_acc', 'State')

  return(first_pass)
}


### PDS_summary

# for each PDS can calculate:
#   1. Rank Abundance
#   2. change from last year
#   3. increasing/decreasing?
#   4. rate of increase?
#   5. Geographic area? (continents?)
#   6. Countries
#   7. Ag Hosts



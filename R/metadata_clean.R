

#' Return a string of names of the matching regex patterns
#'
#' @param pattern_vec a named vector of regex patterns to use
#' @param search_string a sting to search for matches
#'
#' @return names of matching patterns concatenated with '_'
#' @export
#'
#' @examples #soon
return_ag_match <-
  function(pattern_vec, search_string){
    # browser()
  match_vec <- purrr::map(.x =pattern_vec, .f=~grepl(.x, search_string, ignore.case = TRUE)) |>
    unlist()

  res <- names(match_vec)[match_vec] |>
    paste(collapse = '_')
  return(res)
}


#' Extract a consensus ag host species from metadata
#'
#' @param dat an ncbi pathogen detection metadata table
#'
#' @return returns a tibble of 2 columns, 1st = target_acc, 2nd = ag_match
#' @export
#'
#' @examples #soon
extract_consensus_ag_species <- function(dat){
  # browser()
  pattern_vec <-
    c(Swine="swine|pork|porcine|sow|sus|hog|pig|scrofa",
      Bovine="bovine|beef|veal|cow|cattle|bos|steer|taurus|calf|bull|dairy|milk",
      Chicken="chicken|chick|gallus|broiler|egg",
      Turkey="turkey|meleagris|gallopavo",
      Human="human|homo|sapiens",
      Horse="equine|equus|horse|caballus",
      Dog="canine|dog|canis",
      Cat= "\\bcat\\b|felis|catus",
      Goat="caprine|goat",
      Sheep="\\bovine|sheep|lamb",
      Duck="duck",
      Goose='goose')

  first_pass <-
    dat |>
    dplyr::transmute(target_acc,
              search_vals=paste(isolation_source, host, ontological_term,epi_type, sep = '_')) |>
    dplyr::mutate(ag_match=map_chr(.x = search_vals, ~return_ag_match(pattern_vec, search_vec = .x)))

    finished <- first_pass |>
      dplyr::filter(ag_match != '')

    second_pass <- first_pass |>
      dplyr::filter(ag_match == '') |>
      dplyr::mutate(ag_match=ifelse(grepl('clinical', search_vals), 'Human', 'Other'))


    result <- dplyr::bind_rows(finished, second_pass) |>
              dplyr::select(target_acc, ag_match)
    return(result)


}

#
# meta <- read_tsv('~/Documents/O157_overview/data/O157:H7_meta.tsv')
#
# TEST <- extract_consensus_ag_species(meta)
# TEST |> count(ag_match)




#' add a Year column to PDD metadata containing the earliest year from the
#' available 'date' fields.
#'
#' @param PDD_metadata_table
#'
#' @return returns the input metadata table with an added `Year` column
#' @export
#'
#' @examples #soon
get_earliest_year <- function(PDD_metadata_table){
  # given a PDD metadata table, extract the earliest year from the date
  # columns and return the original table with the 'Year' variable added
  result <-
    PDD_metadata_table |>
    dplyr::select(target_acc, ends_with('date')) |>
    dplyr::mutate(across(.cols = ends_with('date'), .fns = as.character)) |>
    tidyr::pivot_longer(cols = ends_with('date'), names_to = 'type', values_to = 'date') |>
    dplyr::mutate(year = as.numeric(sub('([0-9][0-9][0-9][0-9]).*','\\1',date))) |>
    dplyr::group_by(target_acc) |>
    dplyr::summarise(Year=year[which.min(year)]) |>
    dplyr::left_join(PDD_metadata_table)

  return(result)

}



#' Return a string of names of the matching regex patterns
#'
#' @param pattern_vec a named vector of regex patterns to use
#' @param search_string a sting to search for matches
#'
#' @return names of matching patterns concatenated with '_'
#'
#'
#' @examples #pat_vec <- c(Swine='hog|swine|sow'); pdtools:::return_ag_match(pattern_vec=pat_vec, 'Hog')
#'
return_ag_match <-
  function(pattern_vec, search_string){
    # browser()
  match_vec <- purrr::map(.x =pattern_vec, .f=~base::grepl(.x, search_string, ignore.case = TRUE)) |>
    base::unlist()

  res <- names(match_vec)[match_vec] |>
    base::paste(collapse = '_')
  return(res)
}


#' Extract a consensus ag host species from metadata
#'
#' @param dat an ncbi pathogen detection metadata table
#'
#' @return returns a tibble of 2 columns, 1st = target_acc, 2nd = ag_match
#' @export
#'
#' @examples extract_consensus_ag_species(klebsiella_example_dat)
#' @importFrom rlang .data
extract_consensus_ag_species <- function(dat){
  # browser()
  pattern_vec <-
    base::c(Swine="swine|pork|porcine|sow|sus|hog|pig|scrofa",
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
    dplyr::transmute(target_acc=.data$target_acc,
              search_vals=base::paste(.data$isolation_source, .data$host, .data$ontological_term,.data$epi_type, sep = '_')) |>
    dplyr::mutate(ag_match=purrr::map_chr(.x = .data$search_vals, ~return_ag_match(pattern_vec, search_string = .x)))

    finished <- first_pass |>
      dplyr::filter(.data$ag_match != '')

    if (base::nrow(finished) == base::nrow(dat)) {
      result <- finished |> dplyr::select(.data$target_acc, .data$ag_match)
      return(result)
    } else {

      second_pass <- first_pass |>
        dplyr::filter(.data$ag_match == '') |>
        dplyr::mutate(ag_match=base::ifelse(base::grepl('clinical', .data$search_vals), 'Human', 'Other'))


      result <- dplyr::bind_rows(finished, second_pass) |>
        dplyr::select(.data$target_acc, .data$ag_match)

      return(result)

    }




}


#' return a Year column containing the earliest year from the
#' available 'date' fields.
#'
#' @param PDD_metadata_table an ncbi pathogen detection metadata table
#'
#' @return returns the a tibble with a `Year` column,
#' @export
#'
#' @examples return_earliest_year(klebsiella_example_dat)
#' @importFrom rlang .data
return_earliest_year <- function(PDD_metadata_table){

    result <-
    PDD_metadata_table |>
    dplyr::select(.data$target_acc, dplyr::ends_with('date')) |>
    dplyr::mutate(dplyr::across(.cols = dplyr::ends_with('date'), .fns = base::as.character)) |>
    tidyr::pivot_longer(cols = dplyr::ends_with('date'), names_to = 'type', values_to = 'date') |>
    dplyr::mutate(year = base::as.numeric(base::sub('([0-9][0-9][0-9][0-9]).*','\\1',.data$date))) |>
    dplyr::group_by(.data$target_acc) |>
    dplyr::summarise(Year=.data$year[base::which.min(.data$year)])

  return(result)

}

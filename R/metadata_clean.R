

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
  res <- map(.x =pattern_vec, .f=~grepl(.x, search_string, ignore.case = TRUE)) %>%
    unlist() %>%
    names(.)[.] %>%
    paste(collapse = '_')
  return(res)
}


#' Extract a consensus ag host species from metadata
#'
#' @param dat an ncbi pathogen detection metadata table
#'
#' @return # returns a tibble of 2 columns, 1st = target_acc, 2nd = ag_match
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
    dat %>%
    transmute(target_acc,
              search_vals=paste(isolation_source, host, ontological_term,epi_type, sep = '_')) %>%
    mutate(ag_match=map_chr(.x = search_vals, ~return_ag_match(pattern_vec, search_vec = .x)))

    finished <- first_pass %>%
      filter(ag_match != '')

    second_pass <- first_pass %>%
      filter(ag_match == '') %>%
      mutate(ag_match=ifelse(grepl('clinical', search_vals), 'Human', 'Other'))


    result <- bind_rows(finished, second_pass) %>% select(target_acc, ag_match)
    return(result)


}

#
# meta <- read_tsv('~/Documents/O157_overview/data/O157:H7_meta.tsv')
#
# TEST <- extract_consensus_ag_species(meta)
# TEST %>% count(ag_match)

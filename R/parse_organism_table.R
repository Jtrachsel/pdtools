
#' Load the organism summary table
#'
#' Loads a table that lists the organisms in the NCBI pathogens project and the number of isolates
#' associated with each
#'
#' Requires the RSelenium and janitor packages.
#' Rselenium was a bit of a pain for me to get working with firefox in my linux env
#'
#' @param browser 'string designating the browser to use, default is firefox'
#' @param url 'url to extract the table from.  default is https://www.ncbi.nlm.nih.gov/pathogens/organisms/
#'
#' @return a dataframe containing the organism information
#' @export
#'
#' @examples 'get_organisms_table()'
get_organism_table <-
  function(browser='firefox',
           url='https://www.ncbi.nlm.nih.gov/pathogens/organisms/'){
    # Open firefox and extract source
    # this was a real pain to get rselenium working, had to downgrade firefox
    # manually delete some licence files etc.
    if(!rlang::is_installed("janitor") | !rlang::is_installed("RSelenium")){
      rlang::abort(c('required packages not available',
                     'please make sure both the RSelenium and janitor packages are installed'))
    }
    rD <- RSelenium::rsDriver(browser = "firefox", verbose = TRUE)
    remDr <- rD[["client"]]
    remDr$navigate(url)
    Sys.sleep(5) # need to give the page time to load?
    html <- remDr$getPageSource()[[1]]


    # Extract table from source
    org_df <- rvest::read_html(html) %>%
      rvest::html_nodes("table") %>% `[[`(1) %>%
      rvest::html_table() %>%
      janitor::clean_names() %>%
      dplyr::filter(.data$species !='Total') %>%
      dplyr::mutate(dplyr::across(.data$clusters:.data$total_isolates,  ~as.numeric(sub(',','',.x))))

    # Close connection
    remDr$close()
    rD[["server"]]$stop()
    return(org_df)

}


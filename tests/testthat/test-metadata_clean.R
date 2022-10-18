

test_that('matches_from_vector_of_patterns returns expected species',{
  ag_patterns <-
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
  expect_equal(pdtools:::matches_from_vector_of_patterns(pattern_vec = ag_patterns, search_string = 'beef_pork'),
               'Swine_Bovine'
)
  })



test_that('extract_consensus_ag_species returns a tibble',{
  expect_equal(
    dim(extract_consensus_ag_species(klebsiella_example_dat)),
    c(200,2))
  })


test_that('return_earliest_year warns about NA coersions', {
  expect_warning(year_dat <- extract_earliest_year(klebsiella_example_dat))

})

test_that('extract_earliest_year adds an appropriate year column',{
  tst <-function(dat){
    extract_earliest_year(dat) %>% dplyr::pull('Year') %>% is.numeric() %>% unique()
  }
  expect_warning(res <- tst(klebsiella_example_dat))
  expect_true(res)


  })


test_that('extract_collection_agency() extracts the appropriate values',{
  test <- klebsiella_example_dat %>%
          extract_collection_agency() %>%
          dplyr::pull() %>%
          base::unique()
  expect_equal(test, c('CDC'))


})


### furrr functions:
test_that('furrr extract_consensus_ag return the same as the sequential version',{
  future::plan('multisession', workers=2)
  furrr_ag <- extract_consensus_ag_species(klebsiella_example_dat, parallel = TRUE)
  purrr_ag <- extract_consensus_ag_species(klebsiella_example_dat)
  expect_equal(furrr_ag, purrr_ag)
})


test_that('furrr extract_collection_agency return the same as the sequential version',{
  future::plan('multisession', workers=2)
  furrr_collect <- extract_collection_agency(klebsiella_example_dat, parallel = TRUE)
  purrr_collect <- extract_collection_agency(klebsiella_example_dat)
  expect_equal(furrr_collect, purrr_collect)
})


test_that('furrr extract_country return the same as the sequential version',{
  future::plan('multisession', workers=2)
  furrr_collect <- extract_country(klebsiella_example_dat, parallel = TRUE)
  purrr_collect <- extract_country(klebsiella_example_dat)
  expect_equal(furrr_collect, purrr_collect)
})

test_that('extract_state returns an appropriate tibble',{
  test <- extract_state(klebsiella_example_dat)
  expect_equal(dim(test), c(200,2))
  expect_equal(sapply(test, class), c(target_acc='character', State='character'))
})




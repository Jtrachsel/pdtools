

test_that('return_ag_match returns expected species',{
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
  expect_equal(pdtools:::return_ag_match(pattern_vec = ag_patterns, search_string = 'beef_pork'),
               'Swine_Bovine'
)
  })



test_that('extract_consensus_ag_species returns a tibble',{
  expect_equal(
    dim(extract_consensus_ag_species(klebsiella_example_dat)),
    c(200,2))
  })


test_that('return_earliest_year warns about NA coersions', {
  expect_warning(year_dat <- return_earliest_year(klebsiella_example_dat))

})

test_that('return_earliest_year adds an appropriate year column',{
  tst <-function(dat){
    return_earliest_year(dat) |> dplyr::pull(Year) |> is.numeric()
  }
  quiet_tst <- purrr::quietly(tst)
  res <- quiet_tst(klebsiella_example_dat)
  expect_true(res$result)
  expect_true(grepl('NAs introduced', res$warnings))


  })


test_that('extract_collection_agency() extracts the appropriate values',{
  test <- klebsiella_example_dat |>
          extract_collection_agency() |>
          dplyr::pull() |>
          base::unique()
  expect_equal(test, 'CDC')


})

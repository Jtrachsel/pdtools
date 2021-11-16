
test_that('list_PDGs returns a 2 column tibble', {
  skip_if_offline()
  skip_on_cran()
  skip_on_bioc()
  expect_equal(
    ncol(list_PDGs('Campylobacter')), 2
  )
})



test_that('download_most_recent_complete downloads some metadata', {
  skip_if_offline()
  skip_on_cran()
  skip_on_bioc()
  res <- download_most_recent_complete('Serratia')
  dl_files <- list.files(pattern = 'PDG')
  expect_equal(length(dl_files), 2)
  expect_true(file.remove(dl_files[1]))
  expect_true(file.remove(dl_files[2]))

})


test_that('find_most_recent_complete returns a vector of length 2', {
  skip_if_offline()
  skip_on_cran()
  skip_on_bioc()
  expect_equal(
    length(find_most_recent_complete('Serratia')), 2)
})


test_that('download_PDD_metadata downloads some metadata', {
  skip_if_offline()
  skip_on_cran()
  skip_on_bioc()
  PDG <- find_most_recent_complete('Serratia')[1]
  download_PDD_metadata('Serratia', PDG[1])
  dl_files <- list.files(pattern = 'PDG')
  expect_equal(length(dl_files), 2)
  expect_true(file.remove(dl_files[1]))
  expect_true(file.remove(dl_files[2]))

})


test_that('check_complete_PDG returns TRUE for a complete PDG',{
  skip_if_offline()
  skip_on_cran()
  skip_on_bioc()
  # this is circular, the find_most_recent_complete() uses check_complete_PDG...
  PDG <- find_most_recent_complete('Serratia')[1]
  expect_true(check_complete_PDG('Serratia', PDG[1]))


})

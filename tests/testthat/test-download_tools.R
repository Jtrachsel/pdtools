
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



test_that('all snp tree urls match expected pattern', {
  snp_tree_urls <- make_SNPtree_urls(organism = 'Klebsiella', data = klebsiella_example_dat, PDG = 'PDG000000012.1053')
  expect_true(all(grepl('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/PDG[0-9]+\\.[0-9]+/SNP_trees/PDS[0-9]+\\.[0-9]+\\.tar\\.gz', snp_tree_urls)))



})



test_that('list_organisms returns a two column tibble', {
  skip_if_offline()
  orgs <- list_organisms()
  expect_equal(ncol(orgs), 2)
})


test_that('make_ftp_paths returns expected values', {

  test <- make_ftp_paths(klebsiella_example_dat$asm_acc, assembly_summary_path = './kleb_assembly_summary.txt') |>
    dplyr::transmute(asm_acc,
              ftp_path2=ftp_path) |>
    dplyr::left_join(klebsiella_example_dat)

  expect_true(all(test$ftp_path == test$ftp_path2))

})


#
# test_that('make_download_urls returns expected values',{
#   test <- make_download_urls(asm_acc = klebsiella_example_dat$asm_acc,
#                              ftp_paths = klebsiella_example_dat$ftp_path, type='fna')
#
# })

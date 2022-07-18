
# test_that('list_PDGs returns a 2 column tibble', {
#   # skip_if_offline()
#   # skip_on_cran()
#   # skip_on_bioc()
#   expect_equal(
#     ncol(list_PDGs('Campylobacter')), 2
#   )
# })



test_that('download_most_recent_complete downloads some metadata', {
  # skip_if_offline()
  # skip_on_cran()
  # skip_on_bioc()
  res <- download_most_recent_complete('Serratia')
  dl_files <- list.files(pattern = 'PDG')
  expect_equal(length(dl_files), 2)
  expect_true(file.remove(dl_files[1]))
  expect_true(file.remove(dl_files[2]))

})


test_that('find_most_recent_complete returns a vector of length 2', {
  # skip_if_offline()
  # skip_on_cran()
  # skip_on_bioc()
  expect_equal(
    length(find_most_recent_complete('Serratia')), 2)
})


# test_that('download_PDD_metadata downloads some metadata', {
#   skip_if_offline()
#   skip_on_cran()
#   skip_on_bioc()
#   PDG <- find_most_recent_complete('Serratia')[1]
#   download_PDD_metadata('Serratia', PDG[1])
#   dl_files <- list.files(pattern = 'PDG')
#   expect_equal(length(dl_files), 2)
#   expect_true(file.remove(dl_files[1]))
#   expect_true(file.remove(dl_files[2]))
#
# })


test_that('check_complete_PDG returns TRUE for a complete PDG',{
  # skip_if_offline()
  # skip_on_cran()
  # skip_on_bioc()
  # this is circular, the find_most_recent_complete() uses check_complete_PDG...
  PDG <- find_most_recent_complete('Serratia')[1]
  expect_true(check_complete_PDG('Serratia', PDG[1]))


})



test_that('all snp tree urls match expected pattern', {
  snp_tree_urls <- make_SNPtree_urls(organism = 'Klebsiella', data = klebsiella_example_dat, PDG = 'PDG000000012.1053')
  expect_true(all(grepl('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/PDG[0-9]+\\.[0-9]+/SNP_trees/PDS[0-9]+\\.[0-9]+\\.tar\\.gz', snp_tree_urls)))



})


#
# test_that('list_organisms returns a two column tibble', {
#   # skip_if_offline()
#   orgs <- list_organisms()
#   expect_equal(ncol(orgs), 2)
# })


test_that('list_organisms returns an appropriate tibble', {
  test <- list_organisms()
  col_classes <- sapply(test, class) %>% unlist()
  expect_true(nrow(test) > 1)
  expect_true(col_classes[1] == 'character')
  expect_true(any(grepl('POSIX',col_classes)))
})


test_that('list_PDGs returns and appropriate tibble',{
  test <- list_PDGs(organism = 'Salmonella')
  col_classes <- sapply(test, class) %>% unlist()
  expect_true(nrow(test) > 1)
  expect_true(col_classes[1] == 'character')
  expect_true(any(grepl('POSIX',col_classes)))
})




test_that('make_ftp_paths returns expected values', {

  test <- make_ftp_paths(klebsiella_example_dat, assembly_summary_path = './kleb_assembly_summary.txt') |>
    dplyr::transmute(asm_acc,
              ftp_path2=ftp_path) |>
    dplyr::left_join(klebsiella_example_dat)

  expect_true(all(test$ftp_path == test$ftp_path2))

})


# genome download related
test_that('make_download_urls returns correct columns',{
  test <-
    klebsiella_example_dat %>%
    make_download_urls(type = 'fna') %>%
    make_download_urls(type = 'gff') %>%
    make_download_urls(type = 'gbff') %>%
    dplyr::select(ends_with('download'))
  expect_equal(dim(test), c(200,3))
  expect_equal(colnames(test), c('fna_download', 'gff_download', 'gbff_download'))
  })





# genome download related

test_that('make_dest_paths returns a correct columns' ,{
  test <-
    klebsiella_example_dat %>%
    make_dest_paths(type='fna', dest_dir = './') %>%
    dplyr::select(contains('fna'))
  expect_equal(dim(test), c(200,1))
  expect_equal(colnames(test), 'fna_dest')



  })

test_that('supported_download_types returns an error on unsupported type', {
  expect_error(pdtools:::supported_download_types('bad_type'),
               regexp = NULL)
})


test_that('supported_download_types returns all supported types', {
  test <-
    list(supported_download_types('fna'),
       supported_download_types('gbff'),
       supported_download_types('gff'),
       supported_download_types('gtf'),
       supported_download_types('faa'),
       supported_download_types('cds')) %>%
    unlist() %>%
    unique()
  expect_equal(test, c('_genomic.fna.gz', '_genomic.gbff.gz',
                       '_genomic.gff.gz', '_genomic.gtf.gz',
                       '_protein.faa.gz', '_cds_from_genomic.fna.gz'))
})

test_that('check_if_files_exist returns correct columns',{
  test <-
    klebsiella_example_dat %>%
    make_dest_paths(dest_dir = './', type = 'fna') %>%
    check_if_files_exist(type='fna') %>%
    dplyr::select(contains('exist'))
  expect_equal(dim(test), c(200,1))
  expect_false(test[,1] %>% unlist() %>% unique())

})

# test_that('download_gbk_assembly_summary organism returns a correctly formatted
#           url',{
#             expect_error(
#               download_gbk_assembly_summary('TEST', organism = 'TEST'),
#               regexp = '.*cannot open URL.*'
#               )
#           })

#
# test_that('make_download_urls returns expected values',{
#   test <- make_download_urls(asm_acc = klebsiella_example_dat$asm_acc,
#                              ftp_paths = klebsiella_example_dat$ftp_path, type='fna')
#
# })



test_that('generate_genome_vector("genome_1", 20) returns an appropriate tibble',{
  genome <- generate_genome_vector('genome_1', 20)
  expect_equal(dim(genome), c(20,3))
})

test_that('generate_pangenome() returns an appropriate matrix',{
  expect_equal(dim(generate_pangenome(num_genomes = 3, num_genes = 20)), c(20,3))
})


test_that('build_ppanggolin_file_fastas() returns an appropriate tibble',{
  expect_equal(dim(build_ppanggolin_file_fastas(incomplete_genome_paths = './test.fa')), c(1,3))
})

test_that('pan_mat_to_gene_vec_tibble() returns an appropriate tibble',{
  expect_equal(dim(pan_mat_to_gene_vec_tibble(generate_pangenome())), c(100,2))
})

test_that('get_pangenome_representatives() returns an appropriate list',{
  test <- get_pangenome_representatives(pan_mat = generate_pangenome(),
                                        desired_coverage = 1,
                                        verbose = TRUE)
  expect_equal(lapply(test, typeof) %>% unlist(), c('character', 'integer', 'double'))
})

#
# test_that('get_pangenome_representatives2() returns an appropriate list',{
#   gvt <- pan_mat_to_gene_vec_tibble(generate_pangenome())
#   test <- get_pangenome_representatives2(gvt, desired_coverage = 1)
#   expect_equal(lapply(test, typeof) %>% unlist(), c('character', 'integer', 'double'))
# })

test_that('remove_strict_core returns a PA matrix', {
  test <- generate_pangenome(core_genome_fraction = 1)
  core_rem <- test %>% remove_strict_core()
  expect_true(nrow(test) > nrow(core_rem))
})


test_that('mark_outliers returns a correct tibble', {
  pan_PA <- pdtools:::generate_pangenome(core_genome_fraction = 1)
  outliers <-
    pan_PA %>%
    t() %>%
    dist(method = 'binary') %>%
    mark_outliers()
  expect_true(nrow(outliers) == ncol(pan_PA))
})

test_that('cluster_genomes returns an appropriate tibble', {
  pan_PA <- pdtools:::generate_pangenome(core_genome_fraction = 1) %>%
    t()
  test <- cluster_genomes(dat_mat = pan_PA)
  expect_true(nrow(test) == nrow(pan_PA))
})

test_that('pick_derep_sets returns an appropriate tibble',{
  pan_PA <- pdtools:::generate_pangenome(core_genome_fraction = 1)
  tests <- pick_derep_sets(pan_PA, num_sets = 5)
  expect_true(nrow(tests) == 5)

  result_type <- lapply(tests$selection_set, typeof) %>% unlist() %>% unique()
  expect_true(result_type == 'list')
})


test_that('calculate novelty returns an appropriate tibble',{
  pan_PA <- pdtools:::generate_pangenome(core_genome_fraction = 1)
  derep_sets <- pick_derep_sets(pan_PA,desired_coverage = 1, num_sets = 5)
  test <- calculate_novelty(derep_sets)
  expect_true(ncol(test) == 8)
  expect_true(nrow(test) > 0)
})

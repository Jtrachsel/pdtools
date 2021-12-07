
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
  test <- get_pangenome_representatives(pan_mat = generate_pangenome(), desired_coverage = 1)
  expect_equal(lapply(test, typeof) |> unlist(), c('character', 'integer', 'double'))
})


test_that('get_pangenome_representatives2() returns an appropriate list',{
  gvt <- pan_mat_to_gene_vec_tibble(generate_pangenome())
  test <- get_pangenome_representatives2(gvt, desired_coverage = 1)
  expect_equal(lapply(test, typeof) |> unlist(), c('character', 'integer', 'double'))
})

test_that('remove_strict_core returns a PA matrix', {
  test <- generate_pangenome(core_genome_fraction = 1)
  core_rem <- test |> remove_strict_core()
  expect_true(nrow(test) > nrow(core_rem))
})

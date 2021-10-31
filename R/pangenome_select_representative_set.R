
generate_genome_vector <- function(genome_name, num_genes, core_genome_fraction=.75){

  core_genome_size= floor(num_genes * core_genome_fraction)
  accessory_genome_size= num_genes - core_genome_size

  core_genome <- sample(x=c(0,1),
                        size = core_genome_size,
                        replace = T, prob = c(.05, .95))
  accessory_genome <- sample(x=c(0,1),
                             size = accessory_genome_size,
                             replace = T, prob = c(.75,.25))

  gene_names <- paste0('gene_', 1:num_genes)
  gene_presence_values <- c(core_genome, accessory_genome)
  names(gene_presence_values) <- gene_names

  result <-
    tibble(genome_name=genome_name,
           gene_name=names(gene_presence_values),
           gene_presence=gene_presence_values)
  return(result)



}
# generate_genome_vector(genome_name = 'genome_1', num_genes = 5000)

generate_pangenome <- function(num_genomes=100, num_genes=1000, core_genome_fraction=.75){
  genomes <- paste0('genome_', 1:num_genomes)
  pangenome_matrix <-
    map(genomes,
        ~generate_genome_vector(.x, num_genes = num_genes, core_genome_fraction=core_genome_fraction)) %>%
    bind_rows() %>%
    pivot_wider(names_from = gene_name, values_from = gene_presence) %>%
    column_to_rownames(var='genome_name') %>%
    as.matrix()
  return(pangenome_matrix)
}
# generate_pangenome(core_genome_fraction = .75)


# example_dat <-
#   map(genomes, ~generate_genome_vector(.x, num_genes = 1000)) %>%
#   bind_rows() %>%
#   pivot_wider(names_from = gene_name, values_from = gene_presence)
#
# example_mat <- example_dat %>% column_to_rownames(var='genome_name') %>% as.matrix()
#
# hist(colSums(example_mat), breaks = 100)
# hist(rowSums(example_mat), breaks = 100)
#
#
# tot_genomes <- nrow(example_mat)
# desired_coverage=.95
# desired_score <- best_score * desired_coverage

return_set_score <- function(pan_mat, set_size){
  #subset a pangenome to a random collection of a defined size
  # return a score that describes the proportion of pangenomes total gene content
  # contained within the reduced set.

  pan_mat <- pan_mat[,colSums(pan_mat) > 0]
  tot_genomes <- nrow(pan_mat)
  set_indicies <- sample(x = 1:tot_genomes, size = set_size)
  set_mat <- pan_mat[set_indicies,]
  score <- sum(colSums(set_mat) > 1)
  return(list(score=score, set_indicies=set_indicies))
}

get_gene_content_reps <- function(pan_mat, desired_coverage=.95, starting_set_size=1, num_iters_per_size=500){
  # hopefully get smallest set of genomes that gives desired coverage of pangenome
  # browser()
  best_score <- ncol(pan_mat)
  tot_genomes <- nrow(pan_mat)
  desired_score <- best_score * desired_coverage

  print(paste(tot_genomes, 'total genomes'))
  print(paste(best_score, '= best possible score'))
  print(paste(desired_score, ' = desired score'))
  desired_score_reached <- FALSE

  set_size <- starting_set_size

  if(starting_set_size < 2){
    print('must use starting set size of 2 or smaller')
    set_size <- 1
  }
  while(!desired_score_reached){
    # browser()
    set_size=set_size + 1
    print(paste0('using set size ', set_size))
    set_size_results <- map(.x = 1:num_iters_per_size,.f = ~return_set_score(pan_mat = pan_mat, set_size = set_size))

    results_frame <-
      enframe(set_size_results) %>%
      mutate(scores=map_dbl(value, 'score'),
             genome_indicies=map(value, 'set_indicies'),
             set_size=set_size)

    passing_results <- results_frame %>% filter(scores >= desired_score)
    if (nrow(passing_results) > 0){
      desired_score_reached <- TRUE
      return(passing_results %>% arrange(desc(scores)))
    }
    print('desired score not reached, increasing set size')
  }
}

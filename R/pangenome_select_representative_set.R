
#' Generate a genome tibble containing presence/absence of genes in a fake pangenome
#'
#' @param genome_name name of the genome to generate
#' @param num_genes number of genes the genome should contain
#' @param core_genome_fraction fraction of the genome that is part of the 'core pangenome'
#'
#' @return 3 column tibble 1) genome_name; 2) gene_name 3) gene_presence
#' @export
#'
#' @examples #soon
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

#' Generate a random synthetic pangenome gene presence absence matrix
#'
#' @param num_genomes number of genomes the matrix will contain
#'
#' @param num_genes number of genes the matrix will contain
#' @param core_genome_fraction fraction of genes that are part of the core genome
#'
#' @return gene presence absence matrix (0/1), rows are genomes, columns are genes
#' @export
#'
#' @examples #soon
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

#' return the coverage score (proportion coverage) of a random set of genomes
#'
#' @param pan_mat input pangenome PA matrix
#' @param set_size number of genomes to randomly select
#'
#' @return returns a named list of length 2, [[1]] = scores, [[2]] genome indicies
#' @export
#'
#' @examples #soon
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

#' look for small sets of pangenomes that produce the desired gene content coverage of a pangenome
#'
#' @param pan_mat input pangenome PA matrix
#' @param desired_coverage proportion of the pangenome to cover (.95)
#' @param starting_set_size initial number of genomes to use as a set size
#' @param num_iters_per_size number of samples to take for each set size setp
#'
#' @return tibble of sets of genomes that meet the desired coverage threshold, will be at least one set, but maybe more.
#' @export
#'
#' @examples
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


#' Build a ppanggolin file from fastas
#'
#' @param complete_genomes_paths
#' @param incomplete_genome_paths
#'
#' @return a tibble satisfying the ppanggolin file requirements with the complete genome contigs indicated as circular
#' @export
#'
#' @examples #soon
build_ppanggolin_file_fastas <-
  function(complete_genome_paths=NULL,
           incomplete_genome_paths=NULL){
    # browser()
    complete_genomes_table <- NULL
    incomplete_genomes_table <- NULL
    if (is.null(complete_genome_paths) & is.null(incomplete_genome_paths)){

      errorCondition('you must specify either complete_genome_paths or incomplete_genome_paths')
    }

    if (!is.null(complete_genome_paths)){
      # browser()
      complete_genomes_table <-
        tibble(paths=complete_genome_paths,
               ID=sub('(.*)\\.f.*a$','\\1',basename(paths)),
               fasta=map(.x = paths, .f=readDNAStringSet),
               c_names=map(.x=fasta, .f=names),
               c_ids=map(c_names,~sub('(\\w+).*', '\\1', .x)),
               contig_names=map_chr(.x=c_ids, .f=~paste(.x, collapse='\t'))) |>
        select(ID, paths, contig_names)
    }

    if (!is.null(incomplete_genome_paths)){
      incomplete_genomes_table <-
        tibble(paths=incomplete_genome_paths,
               ID=sub('(.*)\\.f.*a$','\\1',basename(paths)),
               contig_names='') %>%
        select(ID, paths, contig_names)
    }

    result <- bind_rows(complete_genomes_table, incomplete_genomes_table)

    return(result)


  }



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

#' THIS IS A BAD WAY TO LOOK FOR REPS look for small sets of pangenomes that produce the desired gene content coverage of a pangenome
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
get_gene_content_reps <-
  function(pan_mat,
           desired_coverage=.95,
           starting_set_size=1,
           num_iters_per_size=500,
           set_size_step=1,
           PARALLEL=TRUE){
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
      set_size=set_size + set_size_step
      print(paste0('using set size ', set_size))

      if (PARALLEL){
        set_size_results <- future_map(progress=TRUE, furrr_options(seed=TRUE),
                                       .x = 1:num_iters_per_size,
                                       .f = ~return_set_score(pan_mat = pan_mat, set_size = set_size))

      } else {
        set_size_results <- map(.x = 1:num_iters_per_size,.f = ~return_set_score(pan_mat = pan_mat, set_size = set_size))
      }




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



#' Convert a pangenome PA matrix to a tibble of gene vectors
#'
#' @param pan_mat presence absence matrix, genes are columns, rows are genomes 1/0
#'
#' @return a tibble with gene vectors for each genome
#' @export
#'
#' @examples #soon
pan_mat_to_gene_vec_tibble <- function(pan_mat){
# browser()
  gene_vec_tibble <-
    apply(pan_mat > 0,
          1,
          function(logical_vec, char_vec){char_vec[logical_vec]},
          colnames(pan_mat)) %>%
    enframe(name = 'genome_name',
            value = 'gene_vec')
  return(gene_vec_tibble)


}

# test <- pan_mat_to_gene_vec_tibble(pan_PA)

# present_features <- function(logical_vec, char_vec){
#   char_vec[logical_vec]
#   }






#' Select a minimum set of genomes that best represent the gene content of the pangenome
#'
#'
#' @param pan_mat
#' @param desired_coverage
#' @param SEED
#'
#' @return
#' @export
#'
#' @examples
get_pangenome_representatives <-
  function(pan_mat, desired_coverage=.95, SEED=3){
    # hopefully get smallest set of genomes that gives desired coverage of pangenome
    # browser()
    genomes <- pan_mat_to_gene_vec_tibble(pan_mat)

    # random starting genome
    rep_genome_index <- sample(1:nrow(pan_mat), size = 1)

    # starting pangenome
    cumulative_pan <- genomes$gene_vec[[rep_genome_index]]
    cumulative_genomes <- genomes$genome_name[[rep_genome_index]]

    # remove starting genome from remaining genomes
    genomes <- genomes[-rep_genome_index,]

    # best score = total number of genes in pangenome
    best_score <- ncol(pan_mat)
    tot_genomes <- nrow(pan_mat)
    desired_score <- best_score * desired_coverage

    print(paste(tot_genomes, 'total genomes'))
    print(paste(best_score, '= best possible score'))
    print(paste(desired_score, ' = desired score'))

    score <- length(cumulative_pan)
    scores <- c(score)
    print(paste0('starting score = ', score))
    while (score < desired_score){

      best_addition_genome <-
        genomes %>%
        mutate(num_new=map_int(.x = gene_vec, .f= ~(sum(!(is.element(.x, cumulative_pan)))))) %>%
        arrange(desc(num_new)) %>%
        slice_head(n = 1)

      cumulative_pan <- c(cumulative_pan, best_addition_genome$gene_vec[[1]]) %>% unique()
      cumulative_genomes <- c(cumulative_genomes, best_addition_genome$genome_name[[1]])
      score <- length(cumulative_pan)
      scores <- c(scores, score)
      print(paste0('new score = ', score))
      proportion_coverages <- scores/best_score
      print(paste0('proportion covered = ', score/best_score))
    }
    return(list(cumulative_genomes, scores, proportion_coverages))
  }



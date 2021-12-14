
#' Build a ppanggolin file from fastas
#'
#' @param complete_genome_paths vector of paths to complete genome assemblies (all contigs are circular)
#' @param incomplete_genome_paths vector of paths to incomplete genome assemblies (all contigs not circular)
#'
#' @return a tibble satisfying the ppanggolin file requirements with the complete genome contigs indicated as circular
#' @export
#'
#' @examples #build_ppanggolin_file_fastas(incomplete_genome_paths=c('genome1.fasta', 'genome2.fasta))
#' @importFrom rlang .data
build_ppanggolin_file_fastas <-
  function(complete_genome_paths=NULL,
           incomplete_genome_paths=NULL){
    # browser()
    complete_genomes_table <- NULL
    incomplete_genomes_table <- NULL
    if (base::is.null(complete_genome_paths) & base::is.null(incomplete_genome_paths)){

      base::errorCondition('you must specify either complete_genome_paths or incomplete_genome_paths')
    }

    if (!base::is.null(complete_genome_paths)){
      # browser()
      complete_genomes_table <-
       tibble::tibble(paths=complete_genome_paths,
                 ID=sub('(.*)\\.f.*a$','\\1',base::basename(.data$paths)),
                 fasta=purrr::map(.x = .data$paths, .f=Biostrings::readDNAStringSet),
                 c_names=purrr::map(.x=.data$fasta, .f=base::names),
                 c_ids=purrr::map(.data$c_names,~base::sub('(\\w+).*', '\\1', .x)),
                 contig_names=purrr::map_chr(.x=.data$c_ids, .f=~base::paste(.x, collapse='\t'))) |>
        dplyr::select(.data$ID, .data$paths, .data$contig_names)
    }

    if (!base::is.null(incomplete_genome_paths)){
      incomplete_genomes_table <-
        tibble::tibble(paths=incomplete_genome_paths,
                       ID=base::sub('(.*)\\.f.*a$','\\1',basename(.data$paths)),
                       contig_names='') |>
        dplyr::select(.data$ID, .data$paths, .data$contig_names)
    }

    result <- dplyr::bind_rows(complete_genomes_table, incomplete_genomes_table)

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
#' @examples generate_genome_vector(genome_name='genome_1', num_genes=2000)
generate_genome_vector <- function(genome_name, num_genes, core_genome_fraction=.75){

  core_genome_size= base::floor(num_genes * core_genome_fraction)
  accessory_genome_size= num_genes - core_genome_size

  core_genome <- base::sample(x=base::c(0,1),
                        size = core_genome_size,
                        replace = T, prob = base::c(.05, .95))
  accessory_genome <- base::sample(x=base::c(0,1),
                             size = accessory_genome_size,
                             replace = T, prob = base::c(.75,.25))

  gene_names <- base::paste0('gene_', 1:num_genes)
  gene_presence_values <- base::c(core_genome, accessory_genome)
  base::names(gene_presence_values) <- gene_names

  result <-
    tibble::tibble(genome_name=genome_name,
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
#' @return gene presence absence matrix (0/1), rows are genes, columns are genomes
#' @export
#'
#' @examples generate_pangenome()
#' @importFrom rlang .data
generate_pangenome <- function(num_genomes=100, num_genes=1000, core_genome_fraction=.75){
  genomes <- base::paste0('genome_', 1:num_genomes)
  pangenome_matrix <-
    purrr::map(genomes,
        ~generate_genome_vector(.x, num_genes = num_genes, core_genome_fraction=core_genome_fraction)) |>
    dplyr::bind_rows() |>
    tidyr::pivot_wider(names_from = .data$gene_name, values_from = .data$gene_presence) |>
    tibble::column_to_rownames(var='genome_name') |>
    base::as.matrix() |>
    t()
  return(pangenome_matrix)
}






#' Convert a pangenome PA matrix to a tibble of gene vectors
#'
#' @param pan_mat presence absence matrix, genes are columns, rows are genomes 1/0
#'
#' @return a tibble with gene vectors for each genome
#' @export
#'
#' @examples #pan_mat_to_gene_vec_tibble(pangenome_presence_absence_matrix)
#'
pan_mat_to_gene_vec_tibble <- function(pan_mat){
  # browser()

  # expects genes as columns...
  gene_vec_tibble <-
    base::apply(pan_mat > 0,
          2,
          function(logical_vec, char_vec){char_vec[logical_vec]},
          base::rownames(pan_mat)) |>
    tibble::enframe(name = 'genome_name',
                    value = 'gene_vec')
  return(gene_vec_tibble)


}



#' Select a minimum set of genomes that best represent the gene content of the pangenome
#'
#'
#' @param pan_mat a presence absence matrix of 1/0, rows are genomes, columns are genes
#' @param desired_coverage proportion of the pangenome's gene content you want the reduced set to contain (.95)
#' @param SEED random seed to use when selecting the first genome of the collection.
#' @param verbose T/F provides updates via print statements
#'
#' @return returns a list of length 3. 1:names of the genomes, 2:scores for each iteration , 3:proportion coverage for each iteration
#' @export
#'
#' @examples #gen_pangenome_representatives(pan_mat)
#' @importFrom rlang .data
get_pangenome_representatives <-
  function(pan_mat, desired_coverage=.95, SEED=3, verbose=FALSE){
    # hopefully get smallest set of genomes that gives desired coverage of pangenome
    # browser()
    genomes <- pan_mat_to_gene_vec_tibble(pan_mat)

    # random starting genome
    rep_genome_index <- base::sample(1:nrow(genomes), size = 1)

    # starting pangenome
    cumulative_pan <- genomes$gene_vec[[rep_genome_index]]
    cumulative_genomes <- genomes$genome_name[[rep_genome_index]]

    # remove starting genome from remaining genomes
    genomes <- genomes[-rep_genome_index,]

    # best score = total number of genes in pangenome
    best_score <- base::nrow(pan_mat)
    tot_genomes <- base::ncol(pan_mat)
    desired_score <- best_score * desired_coverage

    score <- base::length(cumulative_pan)
    scores <- base::c(score)

    if (verbose) {
      print(base::paste(tot_genomes, 'total genomes'))
      print(base::paste(best_score, '= best possible score'))
      print(base::paste(desired_score, ' = desired score'))
      print(base::paste0('starting score = ', score))
    }

    while (score < desired_score){

      # calculates the number of new genes each genome would contribute to the cumulative pangenome

      genomes <-
        genomes |>
        dplyr::mutate(num_new=purrr::map_int(.x = .data$gene_vec, .f= ~(base::sum(!(base::is.element(.x, cumulative_pan)))))) |>
        dplyr::filter(.data$num_new > 0) # removes genomes that do not contribute new information

      # filters the genomes to only those that contain the max number of new genes for that iteration
      # selects a random genome from those that contribute the max number of new genes
      best_addition_genome <-
        genomes |>
        dplyr::filter(.data$num_new == max(.data$num_new)) |>
        dplyr::slice_sample(n = 1)

      cumulative_pan <- base::c(cumulative_pan, best_addition_genome$gene_vec[[1]]) |> base::unique()
      cumulative_genomes <- base::c(cumulative_genomes, best_addition_genome$genome_name[[1]])
      score <- base::length(cumulative_pan)
      scores <- base::c(scores, score)
      proportion_coverages <- scores/best_score

      if (verbose){

        base::print(base::paste0('new score = ', score))
        base::print(base::paste0('proportion covered = ', score/best_score))

      }
    }
    return(base::list(cumulative_genomes, scores, proportion_coverages))
  }

# genomes <- pan_mat_to_gene_vec_tibble(pan_mat)


#' Removes genes present in all genomes from pangenome presence/absence matrix
#'
#' @param pan_PA a pangenome presence absence matrix
#' @param rows_are_genes a logical indicating if genes are rows in the matrix
#'
#' @return returns a pangenome presence/absence matrix with the strict core removed
#' @export
#'
#' @examples generate_pangenome() |> remove_strict_core()
remove_strict_core <- function(pan_PA, rows_are_genes=NULL){
  # check that strict core exists first!


  if (base::is.null(rows_are_genes)){
    base::print('you did not specify if rows or columns are genes')
    rows_are_genes <- base::ifelse(base::nrow(pan_PA) > base::ncol(pan_PA), TRUE, FALSE)
    base::print(base::paste('I think rows_are_genes =', rows_are_genes))
  }

  if (rows_are_genes){
    pan_PA_strict_core_removed <- pan_PA[base::rowSums(pan_PA) != base::ncol(pan_PA),]
  } else {
    pan_PA_strict_core_removed <- pan_PA[,base::colSums(pan_PA) != base::nrow(pan_PA)]
    }
  return(pan_PA_strict_core_removed)

}


#' Get pangenome representatives from a gene_vec_tibble
#'
#' @param gene_vec_tibble an object returned by pan_mat_to_gene_vec_tibble()
#' @param desired_coverage proportion of gene content desired 0-1
#' @param SEED random seed to use (for selecting 1st genome)
#' @param best_possible_score to save time you can pre-calculate the best possible score (total gene content of pangenome)
#'
#' @return a list of 3; list(cumulative_genomes, scores, proportion_coverages)
#' @export
#'
#' @examples generate_pangenome() |> pan_mat_to_gene_vec_tibble() |> get_pangenome_representatives2()
get_pangenome_representatives2 <-
  function(gene_vec_tibble, desired_coverage=.95, SEED=3, best_possible_score=NULL){
    # hopefully get smallest set of genomes that gives desired coverage of pangenome
    # browser()

    # can save time by providing the best possible score
    # which will be total number of genes in not strict core genome
    if(base::is.null(best_possible_score)){
      base::print('calculating total number of genes, you can speed this up if you supply the best_possible_score parameter')
      best_possible_score <-
        purrr::reduce(gene_vec_tibble$gene_vec, ~base::c(.x, .y) |> base::unique()) |>
        base::length()
    }

    # gene_vec_tibble
    # random starting genome
    chosen_genome_index <- base::sample(1:base::nrow(gene_vec_tibble), size = 1)

    # starting pangenome
    cumulative_pan <- gene_vec_tibble$gene_vec[[chosen_genome_index]]
    cumulative_genomes <- gene_vec_tibble$genome_name[[chosen_genome_index]]

    # remove starting genome from remaining genomes
    gene_vec_tibble <- gene_vec_tibble[-chosen_genome_index,]
    # best_score <- gene_vec_tibble$gene_vec |> unique() |> length()
    # best score = total number of genes in pangenome
    best_score <- best_possible_score
    tot_genomes <- base::nrow(gene_vec_tibble)
    desired_score <- best_score * desired_coverage

    base::print(base::paste(tot_genomes, 'total genomes'))
    base::print(base::paste(best_score, '= best possible score'))
    base::print(base::paste(desired_coverage, '= desired coverage'))
    base::print(base::paste(desired_score, '= desired score'))

    score <- base::length(cumulative_pan)
    scores <- base::c(score)
    base::print(base::paste0('starting score = ', score))
    while (score < desired_score){

      # calculates the number of new genes each genome would contribute to the cumulative pangenome
      # filters to only genomes that will contribute the max possible new genes
      # selects a random one (because all that make it through filter will contibute equally).
      best_addition_genome <-
        gene_vec_tibble |>
        dplyr::mutate(num_new=purrr::map_int(.x = .data$gene_vec, .f= ~(base::sum(!(base::is.element(.x, cumulative_pan)))))) |>
        # dplyr::arrange(dplyr::desc(.data$num_new)) |>
        dplyr::filter(.data$num_new == max(.data$num_new)) |>
        dplyr::slice_sample(n = 1)

      # remove selected genome from remaining genomes
      gene_vec_tibble <- gene_vec_tibble |> dplyr::filter(.data$genome_name != best_addition_genome$genome_name)

      cumulative_pan <- base::c(cumulative_pan, best_addition_genome$gene_vec[[1]]) |> base::unique()
      cumulative_genomes <- base::c(cumulative_genomes, best_addition_genome$genome_name[[1]])
      score <- base::length(cumulative_pan)
      scores <- base::c(scores, score)
      base::print(base::paste0('new score = ', score))
      proportion_coverages <- scores/best_score
      print(base::paste0('proportion covered = ', base::round(score/best_score, digits = 3)))
    }
    return(base::list(cumulative_genomes, scores, proportion_coverages))
  }


# artifacts

# return_set_score <- function(pan_mat, set_size){
#   #subset a pangenome to a random collection of a defined size
#   # return a score that describes the proportion of pangenomes total gene content
#   # contained within the reduced set.
#
#   pan_mat <- pan_mat[,colSums(pan_mat) > 0]
#   tot_genomes <- nrow(pan_mat)
#   set_indicies <- sample(x = 1:tot_genomes, size = set_size)
#   set_mat <- pan_mat[set_indicies,]
#   score <- sum(colSums(set_mat) > 1)
#   return(list(score=score, set_indicies=set_indicies))
# }


# THIS IS A BAD WAY TO LOOK FOR REPS look for small sets of pangenomes that produce the desired gene content coverage of a pangenome

# get_gene_content_reps <-
#   function(pan_mat,
#            desired_coverage=.95,
#            starting_set_size=1,
#            num_iters_per_size=500,
#            set_size_step=1,
#            PARALLEL=TRUE){
#     # hopefully get smallest set of genomes that gives desired coverage of pangenome
#     # browser()
#     best_score <- ncol(pan_mat)
#     tot_genomes <- nrow(pan_mat)
#     desired_score <- best_score * desired_coverage
#
#     print(paste(tot_genomes, 'total genomes'))
#     print(paste(best_score, '= best possible score'))
#     print(paste(desired_score, ' = desired score'))
#     desired_score_reached <- FALSE
#
#     set_size <- starting_set_size
#
#     if(starting_set_size < 2){
#       print('must use starting set size of 2 or smaller')
#       set_size <- 1
#     }
#     while(!desired_score_reached){
#       # browser()
#       set_size=set_size + set_size_step
#       print(paste0('using set size ', set_size))
#
#       if (PARALLEL){
#         set_size_results <- future_map(progress=TRUE, furrr_options(seed=TRUE),
#                                        .x = 1:num_iters_per_size,
#                                        .f = ~return_set_score(pan_mat = pan_mat, set_size = set_size))
#
#       } else {
#         set_size_results <- map(.x = 1:num_iters_per_size,.f = ~return_set_score(pan_mat = pan_mat, set_size = set_size))
#       }
#
#
#
#
#       results_frame <-
#         enframe(set_size_results) |>
#         mutate(scores=map_dbl(value, 'score'),
#                genome_indicies=map(value, 'set_indicies'),
#                set_size=set_size)
#
#       passing_results <- results_frame |> filter(scores >= desired_score)
#       if (nrow(passing_results) > 0){
#         desired_score_reached <- TRUE
#         return(passing_results |> arrange(desc(scores)))
#       }
#       print('desired score not reached, increasing set size')
#     }
#   }




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
                 contig_names=purrr::map_chr(.x=.data$c_ids, .f=~base::paste(.x, collapse='\t'))) %>%
        dplyr::select(.data$ID, .data$paths, .data$contig_names)
    }

    if (!base::is.null(incomplete_genome_paths)){
      incomplete_genomes_table <-
        tibble::tibble(paths=incomplete_genome_paths,
                       ID=base::sub('(.*)\\.f.*a$','\\1',basename(.data$paths)),
                       contig_names='') %>%
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
        ~generate_genome_vector(.x, num_genes = num_genes, core_genome_fraction=core_genome_fraction)) %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_wider(names_from = .data$gene_name, values_from = .data$gene_presence) %>%
    tibble::column_to_rownames(var='genome_name') %>%
    base::as.matrix() %>%
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
          base::rownames(pan_mat)) %>%
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
        genomes %>%
        dplyr::mutate(num_new=purrr::map_int(.x = .data$gene_vec, .f= ~(base::sum(!(base::is.element(.x, cumulative_pan)))))) %>%
        dplyr::filter(.data$num_new > 0) # removes genomes that do not contribute new information

      # filters the genomes to only those that contain the max number of new genes for that iteration
      # selects a random genome from those that contribute the max number of new genes
      best_addition_genome <-
        genomes %>%
        dplyr::filter(.data$num_new == max(.data$num_new)) %>%
        dplyr::slice_sample(n = 1)

      cumulative_pan <- base::c(cumulative_pan, best_addition_genome$gene_vec[[1]]) %>% base::unique()
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
#' @examples generate_pangenome() %>% remove_strict_core()
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
#' @examples generate_pangenome() %>% pan_mat_to_gene_vec_tibble() %>% get_pangenome_representatives2()
get_pangenome_representatives2 <-
  function(gene_vec_tibble, desired_coverage=.95, SEED=3, best_possible_score=NULL){
    # hopefully get smallest set of genomes that gives desired coverage of pangenome
    # browser()

    # can save time by providing the best possible score
    # which will be total number of genes in not strict core genome
    if(base::is.null(best_possible_score)){
      base::print('calculating total number of genes, you can speed this up if you supply the best_possible_score parameter')
      best_possible_score <-
        purrr::reduce(gene_vec_tibble$gene_vec, ~base::c(.x, .y) %>% base::unique()) %>%
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
    # best_score <- gene_vec_tibble$gene_vec %>% unique() %>% length()
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
        gene_vec_tibble %>%
        dplyr::mutate(num_new=purrr::map_int(.x = .data$gene_vec, .f= ~(base::sum(!(base::is.element(.x, cumulative_pan)))))) %>%
        # dplyr::arrange(dplyr::desc(.data$num_new)) %>%
        dplyr::filter(.data$num_new == max(.data$num_new)) %>%
        dplyr::slice_sample(n = 1)

      # remove selected genome from remaining genomes
      gene_vec_tibble <- gene_vec_tibble %>% dplyr::filter(.data$genome_name != best_addition_genome$genome_name)

      cumulative_pan <- base::c(cumulative_pan, best_addition_genome$gene_vec[[1]]) %>% base::unique()
      cumulative_genomes <- base::c(cumulative_genomes, best_addition_genome$genome_name[[1]])
      score <- base::length(cumulative_pan)
      scores <- base::c(scores, score)
      base::print(base::paste0('new score = ', score))
      proportion_coverages <- scores/best_score
      print(base::paste0('proportion covered = ', base::round(score/best_score, digits = 3)))
    }
    return(base::list(cumulative_genomes, scores, proportion_coverages))
  }



#' Title
#'
#' @param DIST
#' @param outlier_prob
#'
#' @return
#' @export
#'
#' @examples
mark_outliers <- function(DIST, outlier_prob=.99){

  mainmat <- as.matrix(DIST)
  dists_to_mine <- rowSums(mainmat)
  outliers <- dists_to_mine > quantile(dists_to_mine, probs = outlier_prob)
  print(paste0('detected ',sum(outliers), ' outliers', ' at ', outlier_prob ,' prob'))
  tibble(asm_acc=names(outliers),
         is_outlier=outliers)

}

#' Cluster genomes at 4 levels from a pangenome gene PA matrix
#' Needs parallelDist, igraph,
#'  genomes as rows and genes as columns?
#'
#' @param dat_mat
#' @param scut
#' @param tcut
#' @param qcut
#' @param pcut
#' @param DIST_METHOD
#' @param output_directory
#' @param write_dist
#'
#' @return
#' @export
#'
#' @examples
cluster_genomes <-
  function(dat_mat,
           pcut=0,
           scut=.5,
           tcut=.75,
           qcut=.75,
           DIST_METHOD='simpson',
           output_directory=NULL,
           write_dist=TRUE){
    # browser()
    dist_filename <- paste0(output_directory, DIST_METHOD, '_dist.rds')
    graph_file_name <- paste0(output_directory, DIST_METHOD, '_graph.rds')

    if (!file.exists(graph_file_name)){
      print('calculating distances')

      # write dist objs for visualization
      # %>%


      if(!file.exists(dist_filename)){
        DIST <-  parallelDist::parallelDist(dat_mat, method=DIST_METHOD)
      } else {
        DIST <- read_rds(dist_filename)
      }

      if(write_dist){

        write_rds(DIST, dist_filename)
      }
      # write_rds('./gifrop_out/islands_overlap_dist.rds')

      print('converting to similarities')

      sim_mat <-
        (1 - DIST) %>%
        as.matrix() #%>%
      # Matrix::Matrix(sparse = T)
      rm(DIST)

      print('building graph')
      # maybe errors when no edges have weight 0?
      g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)

      print('writing graph')
      write_rds(g, graph_file_name)
      # igraph::write.graph(g, file = graph_file_name, format = 'edgelist')

    }  else {
      print(paste0('graph already exists, reading graph from', graph_file_name))
      g <- read_rds(graph_file_name)
    }


    if ( any(E(g)[weight < pcut]) ){
      g <- delete_edges(g, E(g)[weight < pcut])
    }

    # any connection = same cluster
    clust1 <- cluster_louvain(g)


    if ( any(E(g)[weight<scut]) ){
      g <- delete_edges(g, E(g)[weight<scut])
    }

    clust2 <- cluster_louvain(g)

    # print('pruning graph, removing edges with overlap coef of less than XXX')

    if ( any(E(g)[weight<tcut]) ){
      g <- delete_edges(g, E(g)[weight<tcut])
    }


    clust3 <- cluster_louvain(g)

    if ( any(E(g)[weight<qcut]) ){
      g <- delete_edges(g, E(g)[weight<qcut])
    }

    clust4 <- cluster_louvain(g)

    clust_info <- tibble(asm_acc = names(membership(clust1)),
                         primary_cluster = membership(clust1),
                         secondary_cluster = membership(clust2),
                         tertiary_cluster = membership(clust3),
                         quat_cluster = membership(clust4))

    return(clust_info)
  }


#' Title
#'
#' @param VECTOR
#'
#' @return
#' @export
#'
#' @examples
selection_orders <-
  function(VECTOR){
    tibble(genome_name=VECTOR,
           selected_order=1:length(VECTOR))
  }


# wrapper function to pick de-replication sets for specified NARMS serotype designations
#' Pick multiple de-replication sets from a pangenome
#'
#' @param pan_PA A pangenome gene presence/absence matrix, genomes as columns,
#' genes as rows
#' @param output_file a path to save the resulting R object to (RDS format)
#' @param num_sets The number of sets to select (25)
#' @param desired_coverage the proportion of genes in the pangenome to cover (.99)
#'
#' @return a tibble with two columns:
#'  1) the random seeds used,
#'  2) list column containing the dereplication sets
#' @export
#'
#' @examples # soon
pick_derep_sets <-
  function(pan_PA, output_file, num_sets=25, desired_coverage=.99){

    if (!file.exists(output_file)){
      # TIC <- tic()

        selection_PA <- pan_PA
        selection_PA <- selection_PA[rowSums(selection_PA) > 0,]

      derep_sets <-
        tibble(seed=seq(1:num_sets),
               # set90=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .90)),
               # set95=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .95)),
               selection_set=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = pan_PA, SEED = .x, desired_coverage = desired_coverage)),)

      saveRDS(derep_sets, output_file)
      print(group)
      # TOC <- toc()
      # print(TOC)


    } else {
      print('specified output file aready exists...returning it')
      derep_sets <- read_rds(output_file)

    }

    return(derep_sets)


  }



#' Calculate genome novelty from a list of selection sets
#'
#' @param selection_set_results
#'
#' @return a dataframe with the novelty scores of all genomes in the list
#' @export
#'
#' @examples # soon
calculate_novelty <-
  function(selection_set_results){
    # browser()
    genome_summary <-
      selection_set_results %>%
      mutate(genome_vectors=map(selection_set, 1),
             selection_orders=map(genome_vectors, selection_orders)) %>%
      select(selection_orders) %>%
      unnest(selection_orders) %>%
      filter(selected_order !=1) %>%
      group_by(genome_name) %>%
      summarise(mean_rank=mean(selected_order),
                median_rank=median(selected_order),
                number_selections=n(),
                best_rank=min(selected_order),
                worst_rank=max(selected_order)) %>%
      arrange(median_rank) %>%
      mutate(novelty_score1=number_selections*((1/(mean_rank + median_rank))),
             novelty_score2=(number_selections + (number_selections/(mean_rank + median_rank))),
             novelty_score3=((number_selections^2 / (mean_rank + median_rank))),
             novelty_score4=((number_selections / (median_rank))),
             novelty_score5=((number_selections^2 / (median_rank)^2)),
             novelty_score6=number_selections/max(number_selections) / (median_rank/n())) %>%

      filter(!(number_selections == 1 & best_rank == 1)) %>%
      transmute(asm_acc=genome_name,
                median_rank,
                number_selections,
                best_rank,
                worst_rank,
                novelty_score=novelty_score4,
                log_novelty=log(novelty_score)) %>%
      arrange(desc(novelty_score)) %>%
      ungroup() %>%
      mutate(RANK=1:n())

    return(genome_summary)

  }




#' Generate a dataframe for plotting the progress of selecting a set of genomes
#'
#' @param result
#' @param seed_num
#'
#' @return
#' @export
#'
#' @examples
res_plot_dat <-
  function(result, seed_num){
    scores <- result[[3]]
    # plot(1:length(scores), scores)
    tibble(seed_num,num_genomes=1:length(scores), scores)
  }

#
#
# grep_vector <- function(pattern_vector, search_vector, INVERT=F){
#
#   matches <- unique(grep(paste(pattern_vector,collapse="|"),
#                          search_vector,invert = INVERT, value=TRUE))
#   return(matches)
#
# }
#



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
#         enframe(set_size_results) %>%
#         mutate(scores=map_dbl(value, 'score'),
#                genome_indicies=map(value, 'set_indicies'),
#                set_size=set_size)
#
#       passing_results <- results_frame %>% filter(scores >= desired_score)
#       if (nrow(passing_results) > 0){
#         desired_score_reached <- TRUE
#         return(passing_results %>% arrange(desc(scores)))
#       }
#       print('desired score not reached, increasing set size')
#     }
#   }



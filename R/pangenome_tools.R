
#' Build a ppanggolin file from fastas
#'
#' @param complete_genome_paths vector of paths to complete genome assemblies (all contigs are circular)
#' @param incomplete_genome_paths vector of paths to incomplete genome assemblies (all contigs not circular)
#'
#' @return a tibble satisfying the ppanggolin file requirements with the complete genome contigs indicated as circular
#' @export
#'
#' @examples build_ppanggolin_file_fastas(incomplete_genome_paths=c('./genomes/genome1.fasta', './genomes/genome2.fasta'))
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
                 # fasta=purrr::map(.x = .data$paths, .f=Biostrings::readDNAStringSet),
                 # c_names=purrr::map(.x=.data$fasta, .f=base::names),
                 c_ids=purrr::map(.data$paths,~get_fasta_contig_names(.x)),
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


#' get contig names from a fasta file
#'
#' @param path path to a fasta file
#'
#' @return a vector of contig names
#' @noRd
#'
#'#'@examples
#'\dontrun{
#'get_fasta_contig_names('genome_1.fasta')
#'}
get_fasta_contig_names <- function(path){

  con <- base::file(path, "r")
  lines <- base::c()
  while(TRUE) {
    line = base::readLines(con, 1)
    if(base::length(line) == 0) break
    else if(base::grepl("^>", line)){
      line <- base::sub('>','',line)
      lines <- base::c(lines, line)

      }
  }
  base::close.connection(con)
  return(lines)
}



#' Generate a genome tibble containing presence/absence of genes in a fake pangenome
#'
#' @param genome_name name of the genome to generate
#' @param num_genes number of genes the genome should contain
#' @param core_genome_fraction fraction of the genome that is part of the 'core pangenome'
#'
#' @return 3 column tibble 1) genome_name; 2) gene_name 3) gene_presence
#' @noRd
#'
#' @examples #genome1_vec <- generate_genome_vector(genome_name='genome_1', num_genes=2000)
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

#' Generate a random synthetic pangenome gene presence absence matrix
#'
#' @param num_genomes number of genomes the matrix will contain
#'
#' @param num_genes number of genes the matrix will contain
#' @param core_genome_fraction fraction of genes that are part of the core genome
#'
#' @return gene presence absence matrix (0/1), rows are genes, columns are genomes
#' @noRd
#'
#' @examples
#' pangenome <- generate_pangenome()
#' pangenome[1:5, 1:5]
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
#' @noRd
#'
#' @examples
#' gvt <- pan_mat_to_gene_vec_tibble(example_pangenome_matrix)
#' gvt[1:5,]
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
#' @examples
#' # this example pangenome has 100 genomes with 1000 total genes
#' # ~5 genomes can provide > 95% gene coverage
#' pan_reps <- get_pangenome_representatives(example_pangenome_matrix)
#' pan_reps
#'
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


      if (verbose){

        base::print(base::paste0('new score = ', score))
        base::print(base::paste0('proportion covered = ', score/best_score))

      }
    }
    # could move this outside while?
    proportion_coverages <- scores/best_score

    return(base::list(cumulative_genomes, scores, proportion_coverages))
  }


#' Removes genes present in all genomes from pangenome presence/absence matrix
#'
#' @param pan_PA a pangenome presence absence matrix
#' @param rows_are_genes a logical indicating if genes are rows in the matrix
#'
#' @return returns a pangenome presence/absence matrix with the strict core removed
#' @export
#'
#' @examples
#' dim(example_pangenome_matrix)
#' pan_mat <- example_pangenome_matrix %>% remove_strict_core()
#' dim(pan_mat)
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


#' #' Get pangenome representatives from a gene_vec_tibble
#' #'
#' #' @param gene_vec_tibble an object returned by pan_mat_to_gene_vec_tibble()
#' #' @param desired_coverage proportion of gene content desired 0-1
#' #' @param SEED random seed to use (for selecting 1st genome)
#' #' @param best_possible_score to save time you can pre-calculate the best possible score (total gene content of pangenome)
#' #' @param max_genomes Maximum number of genomes to select
#' #'
#' #' @return a list of 3; list(cumulative_genomes, scores, proportion_coverages)
#' #' @noRd
#' #'
#' #' @examples #generate_pangenome() %>% pan_mat_to_gene_vec_tibble() %>% get_pangenome_representatives2()
#' get_pangenome_representatives2 <-
#'   function(gene_vec_tibble,
#'            desired_coverage=.95,
#'            SEED=3,
#'            best_possible_score=NULL,
#'            max_genomes=1000){
#'     # hopefully get smallest set of genomes that gives desired coverage of pangenome
#'     # browser()
#'
#'     # can save time by providing the best possible score
#'     # which will be total number of genes in not strict core genome
#'     if(base::is.null(best_possible_score)){
#'       base::print('calculating total number of genes, you can speed this up if you supply the best_possible_score parameter')
#'       best_possible_score <-
#'         purrr::reduce(gene_vec_tibble$gene_vec, ~base::c(.x, .y) %>% base::unique()) %>%
#'         base::length()
#'     }
#'
#'     # gene_vec_tibble
#'     # random starting genome
#'     chosen_genome_index <- base::sample(1:base::nrow(gene_vec_tibble), size = 1)
#'
#'     # starting pangenome
#'     cumulative_pan <- gene_vec_tibble$gene_vec[[chosen_genome_index]]
#'     cumulative_genomes <- gene_vec_tibble$genome_name[[chosen_genome_index]]
#'
#'     # remove starting genome from remaining genomes
#'     gene_vec_tibble <- gene_vec_tibble[-chosen_genome_index,]
#'     # best_score <- gene_vec_tibble$gene_vec %>% unique() %>% length()
#'     # best score = total number of genes in pangenome
#'     best_score <- best_possible_score
#'     tot_genomes <- base::nrow(gene_vec_tibble)
#'     desired_score <- best_score * desired_coverage
#'
#'     base::print(base::paste(tot_genomes, 'total genomes'))
#'     base::print(base::paste(best_score, '= best possible score'))
#'     base::print(base::paste(desired_coverage, '= desired coverage'))
#'     base::print(base::paste(desired_score, '= desired score'))
#'
#'     score <- base::length(cumulative_pan)
#'     scores <- base::c(score)
#'     base::print(base::paste0('starting score = ', score))
#'     while (score < desired_score & length(cumulative_genomes) < max_genomes){
#'
#'       # calculates the number of new genes each genome would contribute to the cumulative pangenome
#'       # filters to only genomes that will contribute the max possible new genes
#'       # selects a random one (because all that make it through filter will contibute equally).
#'       best_addition_genome <-
#'         gene_vec_tibble %>%
#'         dplyr::mutate(num_new=purrr::map_int(.x = .data$gene_vec, .f= ~(base::sum(!(base::is.element(.x, cumulative_pan)))))) %>%
#'         # dplyr::arrange(dplyr::desc(.data$num_new)) %>%
#'         dplyr::filter(.data$num_new == max(.data$num_new)) %>%
#'         dplyr::slice_sample(n = 1)
#'
#'       # remove selected genome from remaining genomes
#'       gene_vec_tibble <- gene_vec_tibble %>% dplyr::filter(.data$genome_name != best_addition_genome$genome_name)
#'
#'       cumulative_pan <- base::c(cumulative_pan, best_addition_genome$gene_vec[[1]]) %>% base::unique()
#'       cumulative_genomes <- base::c(cumulative_genomes, best_addition_genome$genome_name[[1]])
#'       score <- base::length(cumulative_pan)
#'       scores <- base::c(scores, score)
#'       base::print(base::paste0('new score = ', score))
#'       proportion_coverages <- scores/best_score
#'       print(base::paste0('proportion covered = ', base::round(score/best_score, digits = 3)))
#'     }
#'     return(base::list(cumulative_genomes, scores, proportion_coverages))
#'   }



#' Mark outliers from a distance matrix
#'
#' @param DIST a dist object
#' @param outlier_prob probability to consider an entity an outlier
#'
#' @return a tibble with an outlier designation for each genome
#' @export
#'
#' @examples
#' mark_outliers(example_pan_dist, outlier_prob=.95)
mark_outliers <- function(DIST, outlier_prob=.99){

  mainmat <- base::as.matrix(DIST)
  dists_to_mine <- base::rowSums(mainmat)
  outliers <- dists_to_mine > stats::quantile(dists_to_mine, probs = outlier_prob)
  base::print(base::paste0('detected ',base::sum(outliers), ' outliers', ' at ', outlier_prob ,' prob'))
  tibble::tibble(asm_acc=names(outliers),
                 is_outlier=outliers)

}

#' Cluster genomes at 4 levels from a pangenome gene PA matrix
#' Needs parallelDist, igraph,
#'  genomes as rows and genes as columns?
#'
#' @param dat_mat pan genome presence absence matrix (rows are genes)
#' @param scut weight for graph pruning at 2nd level
#' @param tcut weight for graph pruning at 3rd level
#' @param qcut weight for graph pruning at 4th level
#' @param pcut weight for graph pruning at 1st level
#' @param DIST_METHOD distance method for building graph (passed to parallelDist())
#' @param output_directory output directory to save the graph in
#' @param write_dist logical, should the dist object be written?
#'
#' @return returns a tibble with a cluster designation for each genome at 4 levels
#' @export
#'
#' @examples
#' # note the t() because the dist() function finds dists between rows of a matrix
#' cluster_genomes(dat_mat=t(example_pangenome_matrix),
#'                 pcut=.80,
#'                 scut=.85,
#'                 tcut=.90,
#'                 qcut=.95,
#'                 write_dist=FALSE,
#'                 write_graph=FALSE)
#'
#'
#'
cluster_genomes <-
  function(dat_mat,
           pcut=0,
           scut=.5,
           tcut=.75,
           qcut=.75,
           DIST_METHOD='binary',
           output_directory=NULL,
           write_dist=TRUE,
           write_graph=TRUE){
    # browser()
    dist_filename <- base::paste0(output_directory, DIST_METHOD, '_dist.rds')
    graph_file_name <- base::paste0(output_directory, DIST_METHOD, '_graph.rds')

    if (!base::file.exists(graph_file_name)){
      base::print('calculating distances')

      # write dist objs for visualization
      # %>%


      if(!base::file.exists(dist_filename)){
        DIST <-  parallelDist::parallelDist(dat_mat, method=DIST_METHOD)
      } else {
        DIST <- readr::read_rds(dist_filename)
      }

      if(write_dist){

        readr::write_rds(DIST, dist_filename)
      }
      # write_rds('./gifrop_out/islands_overlap_dist.rds')

      base::print('converting to similarities')

      sim_mat <-
        (1 - DIST) %>%
        base::as.matrix() #%>%
      # Matrix::Matrix(sparse = T)
      base::rm(DIST)

      base::print('building graph')
      # maybe errors when no edges have weight 0?
      g <- igraph::graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)
      if(write_graph == TRUE){
        base::print('writing graph')
        readr::write_rds(g, graph_file_name)
      }

      # igraph::write.graph(g, file = graph_file_name, format = 'edgelist')

    }  else {
      base::print(paste0('graph already exists, reading graph from', graph_file_name))
      g <- readr::read_rds(graph_file_name)
    }

    # 1st level:
    bad_edges <- E(g)[E(g)$weight < pcut]
    g <- igraph::delete_edges(g, bad_edges)

    clust1 <- igraph::cluster_louvain(g)

    # 2nd level
    bad_edges <- E(g)[E(g)$weight < scut]
    g <- igraph::delete_edges(g, bad_edges)

    clust2 <- igraph::cluster_louvain(g)

    # 3rd level
    bad_edges <- E(g)[E(g)$weight < tcut]
    g <- igraph::delete_edges(g, bad_edges)

    clust3 <- igraph::cluster_louvain(g)

    # 4th level
    bad_edges <- E(g)[E(g)$weight < qcut]
    g <- igraph::delete_edges(g, bad_edges)

    clust4 <- igraph::cluster_louvain(g)



    clust_info <- tibble::tibble(asm_acc = base::names(igraph::membership(clust1)),
                                 primary_cluster = igraph::membership(clust1),
                                 secondary_cluster = igraph::membership(clust2),
                                 tertiary_cluster = igraph::membership(clust3),
                                 quat_cluster = igraph::membership(clust4))

    return(clust_info)
  }


#' helper function to add selection orders to a vector of genomes
#'
#' @param VECTOR a vector of genome names
#'
#' @return a tibble of two columns, 1=genome_name, 2=selected_order
#' @noRd
#'
#' @examples # selection_orders(genome_vector)
selection_orders <-
  function(VECTOR){
    tibble::tibble(genome_name=VECTOR,
                   selected_order=1:length(VECTOR))
  }



#' Pick multiple de-replication sets from a pangenome, uses furrr::future_map
#' Be sure to set your 'plan()'!!!
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
  function(pan_PA, output_file, num_sets=25, desired_coverage=.95){

    if (!base::file.exists(output_file)){
      # TIC <- tic()

        selection_PA <- pan_PA
        selection_PA <- selection_PA[base::rowSums(selection_PA) > 0,]

      derep_sets <-
        tibble::tibble(seed=base::seq(1:num_sets),
               # set90=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .90)),
               # set95=future_map(.x = seed, .options = furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = selection_PA, SEED = .x, desired_coverage = .95)),
               selection_set=furrr::future_map(.x = .data$seed, .options = furrr::furrr_options(seed = 1), ~ get_pangenome_representatives(pan_mat = pan_PA, SEED = .x, desired_coverage = desired_coverage)),)

      base::saveRDS(derep_sets, output_file)

      # TOC <- toc()
      # print(TOC)


    } else {
      base::print('specified output file aready exists...returning it')
      derep_sets <- readr::read_rds(output_file)

    }

    return(derep_sets)


  }



#' Calculate genome 'novelty' from a list of selection sets returned by pick_derep_sets
#'
#' @param selection_set_results a tibble of selection sets
#'
#' @return a dataframe with the novelty scores of all genomes in the list
#'         'Novelty' score for each genome is (number of selections) / median rank of selection
#' @export
#'
#' @examples # soon
calculate_novelty <-
  function(selection_set_results){
    # browser()
    genome_summary <-
      selection_set_results %>%
      dplyr::mutate(genome_vectors=purrr::map(.data$selection_set, 1),
             selection_orders=purrr::map(.data$genome_vectors, selection_orders)) %>%
      dplyr::select(.data$selection_orders) %>%
      tidyr::unnest(.data$selection_orders) %>%
      dplyr::filter(.data$selected_order !=1) %>%
      dplyr::group_by(.data$genome_name) %>%
      dplyr::summarise(mean_rank=mean(.data$selected_order),
                       median_rank=stats::median(.data$selected_order),
                       number_selections=dplyr::n(),
                       best_rank=base::min(.data$selected_order),
                       worst_rank=base::max(.data$selected_order)) %>%
      dplyr::arrange(.data$median_rank) %>%
      dplyr::mutate(novelty_score1=.data$number_selections*((1/(.data$mean_rank + .data$median_rank))),
             novelty_score2=(.data$number_selections + (.data$number_selections/(.data$mean_rank + .data$median_rank))),
             novelty_score3=((.data$number_selections^2 / (.data$mean_rank + .data$median_rank))),
             novelty_score4=((.data$number_selections / (.data$median_rank))),
             novelty_score5=((.data$number_selections^2 / (.data$median_rank)^2)),
             novelty_score6=.data$number_selections/base::max(.data$number_selections) / (.data$median_rank/dplyr::n())) %>%

      dplyr::filter(!(.data$number_selections == 1 & .data$best_rank == 1)) %>%
      dplyr::transmute(asm_acc=.data$genome_name,
                       .data$median_rank,
                       .data$number_selections,
                       .data$best_rank,
                       .data$worst_rank,
                novelty_score=.data$novelty_score4,
                log_novelty=base::log(.data$novelty_score)) %>%
      dplyr::arrange(dplyr::desc(.data$novelty_score)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(RANK=1:dplyr::n())

    return(genome_summary)

  }




#' Generate a dataframe for plotting the progress of selecting a set of genomes
#'
#' @param result TODO
#' @param seed_num TODO
#'
#' @return returns a tibble containing the score of the set with each additional genome
#'
#'
#' @examples # res_plot_dat()
res_plot_dat <-
  function(result, seed_num){
    scores <- result[[3]]
    # plot(1:length(scores), scores)
    tibble::tibble(seed_num,num_genomes=1:base::length(scores), scores)
  }




#' Get a subset of a pangenome that tries to represent the gene presence absence diversity
#' Selects a random genome,
#' then selects the most distant genome from the selected genome
#' then
#'
#' @param pan_mat gene presence absence matrix
#' @param pan_dist distance matrix if not providing the PA matrix
#' @param SEED random seed
#' @param verbose include print statemetns?
#' @param CUTOFF stop choosing new genomes when distances are below this level
#' @param max_genomes stop choosing new genomes when this many genomes are chosen
#'
#' @return returns a tibble, 1) asm_acc, 2) min_jacc, 3) order
#' @export
#'
#' @examples # get_pangenome_representatives_jaccard(pan_dist)
get_pangenome_representatives_jaccard <-
  function(pan_mat=NULL,pan_dist=NULL,
           SEED=3,
           verbose=FALSE,
           CUTOFF=.50, # quantile of min_jaccs select until
           max_genomes=1000){
    # hopefully get smallest set of genomes that gives desired coverage of pangenome
    # browser()
    base::set.seed(SEED)

    # setDTthreads(threads = 1, restore_after_fork = TRUE)

    if(base::is.null(pan_dist)){
      base::print('calculating distances')
      pan_dist <- parallelDist::parallelDist(base::t(pan_mat), method = 'binary')
      base::print('done with distance calculation')

    }

    # pick a starting genome

    rep_genome <- base::sample(attr(pan_dist,"Labels"), size = 1)
    selected_genomes <- tibble::tibble(to=rep_genome,
                               min_jacc=100)

    selected_genome <- selected_genomes

    # move the distance matrix to a long format
    base::print('converting to long')
    test_long <-
      pan_dist %>%
      base::as.matrix() %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var='from') %>%
      tidyr::pivot_longer(cols=-.data$from, names_to = 'to', values_to = 'jacc') %>%
      dplyr::filter(.data$from != .data$to)
    base::print('done converting to long')

    # calculate the cutoff for when to stop adding genomes to the set
    base::print('calculating cutoff')

    # these are nearest neighbors
    min_jaccs <-
      test_long %>%
      dplyr::group_by(.data$to) %>%
      dplyr::summarise(min_jacc=base::min(.data$jacc))

    cutoff <- stats::quantile(min_jaccs$min_jacc, probs = CUTOFF)
    base::print(base::paste0('cutoff is: ',cutoff))
    base::rm(min_jaccs)

    # because we will stop when the min_jacc to a selected genome is less than this cutoff, we can
    # ignore any jacc distances lower than this
    # og <- nrow(test_long)
    # test_long <- test_long %>% filter.(from != to & jacc > cutoff)
    # og - nrow(test_long)

    # TOCS <- list()
    while (base::nrow(selected_genomes) < max_genomes &
           # is.na(selected_genome$min_jacc) | # only for 1st 'selected_genome'
           selected_genome$min_jacc > cutoff){
      # tic()
      base::print(base::paste0('selecting genome ', base::nrow(selected_genomes)))
      # browser()

      # I think uncommenting this line makes this much slower???
      # now that we've selected a genome, we dont have to consider it as an option
      # in future iterations.  Should reduce the # of rows we are manipulating by n each iteration
      # test_long <- test_long %>% filter(!(to %in% selected_genomes$to))

      # what about this: just remove specifically the last genome selected
      # instead of checking every genome in the current selection set?
      #HERE
      # test_long <- test_long %>% filter(to != selected_genome$to)

      selected_genome <-
        test_long %>%
        dplyr::filter(.data$from %in% selected_genomes$to & #HERE
                 !(.data$to %in% selected_genomes$to)
        ) %>%
        dplyr::group_by(.data$to) %>%
        dplyr::summarise(
          # sum_jacc=sum(jacc), # sum of distances from chosen genome to all selected genomes
          min_jacc=base::min(.data$jacc)) %>%
        dplyr::arrange(dplyr::desc(.data$min_jacc)) %>%
        # arrange(desc(!!rlang::sym(SCORE))) %>%
        # mutate(var_in_top5)
        dplyr::slice_head(n=1) #%>%
      # pull(to)

      selected_genomes <- dplyr::bind_rows(selected_genomes, selected_genome)
      # TOCS[[nrow(selected_genomes)]] <- toc()

    }
    selected_genomes <-
      selected_genomes %>%
      dplyr::ungroup() %>%
      dplyr::transmute(asm_acc=.data$to, .data$min_jacc, ORDER=0:(base::sum(dplyr::n()) - 1))

    # return(selected_genomes)
    return(selected_genomes)
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



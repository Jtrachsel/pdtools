
#' Build a ppanggolin file from fastas
#'
#' @param complete_genomes_paths vector of paths to complete genome assemblies (all contigs are circular)
#' @param incomplete_genome_paths vector of paths to incomplete genome assemblies (all contigs not circular)
#'
#' @return a tibble satisfying the ppanggolin file requirements with the complete genome contigs indicated as circular
#' @export
#'
#' @examples build_ppanggolin_file_fastas(complete_genomes=complete_genome_paths, incomplete_genomes=incomplete_genome_paths)
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
       tibble::tibble(paths=complete_genome_paths,
                 ID=sub('(.*)\\.f.*a$','\\1',basename(paths)),
                 fasta=purrr::map(.x = paths, .f=Biostrings::readDNAStringSet),
                 c_names=purrr::map(.x=fasta, .f=names),
                 c_ids=purrr::map(c_names,~sub('(\\w+).*', '\\1', .x)),
                 contig_names=purrr::map_chr(.x=c_ids, .f=~paste(.x, collapse='\t'))) |>
        select(ID, paths, contig_names)
    }

    if (!is.null(incomplete_genome_paths)){
      incomplete_genomes_table <-
        tibble::tibble(paths=incomplete_genome_paths,
                       ID=sub('(.*)\\.f.*a$','\\1',basename(paths)),
                       contig_names='') |>
        dplyr::select(ID, paths, contig_names)
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
#' @return gene presence absence matrix (0/1), rows are genomes, columns are genes
#' @export
#'
#' @examples generate_pangenome()
generate_pangenome <- function(num_genomes=100, num_genes=1000, core_genome_fraction=.75){
  genomes <- paste0('genome_', 1:num_genomes)
  pangenome_matrix <-
    purrr::map(genomes,
        ~generate_genome_vector(.x, num_genes = num_genes, core_genome_fraction=core_genome_fraction)) |>
    dplyr::bind_rows() |>
    tidyr::pivot_wider(names_from = gene_name, values_from = gene_presence) |>
    tibble::column_to_rownames(var='genome_name') |>
    as.matrix()
  return(pangenome_matrix)
}






#' Convert a pangenome PA matrix to a tibble of gene vectors
#'
#' @param pan_mat presence absence matrix, genes are columns, rows are genomes 1/0
#'
#' @return a tibble with gene vectors for each genome
#' @export
#'
#' @examples pan_mat_to_gene_vec_tibble(pangenome_presence_absence_matrix)
pan_mat_to_gene_vec_tibble <- function(pan_mat){
# browser()
  gene_vec_tibble <-
    apply(pan_mat > 0,
          1,
          function(logical_vec, char_vec){char_vec[logical_vec]},
          colnames(pan_mat)) |>
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
#'
#' @return returns a list of length 3. [[1]]=names of the genomes, [[2]]=scores for each iteration , [[3]]=proportion coverage for each iteration
#' @export
#'
#' @examples gen_pangenome_representatives(pan_mat)
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

      # calculates the number of new genes each genome would contribute to the cumulative pangenome
      # arranges the genomes by the number of new genes they would contribute
      # selects the first one and adds it to the cumulative pangenome.
      best_addition_genome <-
        genomes |>
        dplyr::mutate(num_new=purrr::map_int(.x = gene_vec, .f= ~(sum(!(is.element(.x, cumulative_pan)))))) |>
        dplyr::arrange(desc(num_new)) |>
        dplyr::slice_head(n = 1)

      cumulative_pan <- c(cumulative_pan, best_addition_genome$gene_vec[[1]]) |> unique()
      cumulative_genomes <- c(cumulative_genomes, best_addition_genome$genome_name[[1]])
      score <- length(cumulative_pan)
      scores <- c(scores, score)
      print(paste0('new score = ', score))
      proportion_coverages <- scores/best_score
      print(paste0('proportion covered = ', score/best_score))
    }
    return(list(cumulative_genomes, scores, proportion_coverages))
  }



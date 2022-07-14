#' NCBI Pathogen Detection metadata for 200 Klebsiella isolates
#'
#' An example dataset of 100 Klebsiella isolates included for testing and
#' example purposes
#'
#' @format A data frame with 200 rows and 64 variables:
#'
#'
#' @source \url{https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/}
"klebsiella_example_dat"



#' A named vector of country names
#'
#' Used to help extract country of isolation from geo_loc_tag column,
#' so far only US and UK have special entries.  Taken from the maps package and slightly modified.
#'
#'
#' @format A named vector of regex patterns (only US and UK have special entries)
#'
#'
#' @source the maps package
"country_vector"

#' A gene presence absence matrix for a synthetic pangenome
#'
#' Used for pangenome function examples
#'
#'
#' @format A presence/absence matrix with genes as rows and genomes as columns
#'
#'
#' @source pdtools:::generate_pangenome()
"example_pangenome_matrix"

#' A distance matrix describing the jaccard distance between genomes in the
#' provided example pangenome
#'
#' Used for pangenome function examples
#'
#'
#' @format A dist object
#'
#'
#' @source dist(method = 'binary', t(example_pangenome_matrix))

"example_pan_dist"

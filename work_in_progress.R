library(pdtools)
library(tidyverse)
library(rvest)
library(lubridate)
library(RCurl)
library(curl)

PDGs <- list_PDGs('Campylobacter')


# download most recent complete data
download_most_recent_complete('Campylobacter')
# PDG000000003.1527


library(data.table)
camp <- fread('./PDG000000003.1527.amr.metadata.tsv', quote='')
camp$isolation_source %>% unique()

camp <- camp %>% filter(host != 'Homo sapiens')
camp <- camp %>% filter(host != 'Chicken')
camp_ovine <- camp %>% filter(grepl('Sheep|Lamb|Ovine', isolation_source))

camp$host %>% table()

camp$isolation_source %>% unique()


# ideas



#
# camp$isolation_source %>% unique()
# camp$host %>% unique()
# camp$ontological_term
#

# Canine
# Feline
# Human
# Primate_other
# Swine
# Chicken
# Turkey
# Bovine
# SHeep
# Goat
# Horse
# Rabbit
# Rodent#?
# Bird_other
# Bovine
# Reptile



#function to check for number of null values in an isolate?

#function to extract broad host from :
# isolation_source
# host
# ontol




tst <-
  generate_pangenome(core_genome_fraction = .1, num_genomes = 1000, num_genes = 10000) %>%
  get_gene_content_reps(desired_coverage = .95, starting_set_size = 10, num_iters_per_size = 250)





param_sweep <- tibble(core_genome_fraction=seq(from=.01, to =.5, by = .01),
                      pangenome=list(generate_pangenome(core_genome_fraction = core_genome_fraction)),
                      results=map(pangenome, ~get_gene_content_reps(.x)))

param_sweep %>% map(results)


param_sweep$results

tst$value
hist(tst$scores)




# get_gene_content_reps(pan_mat = example_mat, desired_coverage = .99, starting_set_size = 5)

# return_set_score(example_mat, set_size = 10)
# return_set_score(example_mat, set_size = 10)
# return_set_score(example_mat, set_size = 10)
#



#
#
# example_mat <- example_mat[,colSums(example_mat) > 0]
# best_score <- ncol(example_mat)

# if desired coverage is reached, reduce the size of the set and recalc
# in progress..
# scores <- c()
# num_iter <- 100
# for (x in 1:num_iter) {
#   genome_set <- sample(1:tot_genomes, size = 5)
#   set_mat <- example_mat[genome_set,]
#   score <- sum(colSums(set_mat) > 1)
#   scores <- c(scores, score)
#
# }
# return(scores)
#
# hist(scores)
#
#
# now i want to select a set of genomes that maximizes the coverage of the pangenome
# I want to remove redundant genomes, select the smallest set of genomes that represents
# the maximum gene content of the pangenome




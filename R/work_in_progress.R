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

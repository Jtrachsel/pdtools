# #Metadata cleanup functions
# install.packages("ritis")
#
#
# library(ritis)
#
# publications(tsn = 70340)


extract_ag_species <- function(column){

  swinenames <- c('swine','pork','porcine','sow','sus','hog','pig', 'scrofa')
  bovinenames <- c('bovine','beef','veal','cow','cattle','bos','steer', 'taurus', 'calf')
  chickennames <- c('chicken', 'chick', 'gallus', 'broiler', 'egg')
  turkeynames <- c('turkey', 'meleagris', 'gallopavo')
  humannames <- c('human', 'homo', 'sapiens')
  horsenames<- c('equine', 'equus', 'horse', 'caballus')
  dognames<- c('canine', 'dog', 'canis')
  catnames<- c('cat', 'felis', 'catus')
  goatnames <- c('caprine', 'goat','','','')
  sheepnames <- c('ovine', 'sheep', 'lamb')
  ducknames <- c('duck','','')
  goosenames <- c('','','')


}


cleaner <-
  function(column, matches_vector, return_value){
    PATTERN <- paste(matches_vector, collapse = '|')
    MATCHES <- grepl(pattern = PATTERN, x = column)
    column[MATCHES] <- return_value
    return(column)
  }

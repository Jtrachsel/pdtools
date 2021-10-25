# pdtools
Functions for working with NCBI's Pathogen Detection data  

This package is designed to help with downloading a collection of isolates from the NCBI pathogen detection database 


# TODO  
- package dependencies are specified wrong?  
- parse metadata functions  
- update_collection() function?
  - Should take and old metadata file and a new metadata file as inputs.  
  - return a vector of new genomes that were not present in the old list.  
- broad_animal_source() function?  
  - Should look in the metadata for key words and then return standardized ag species:
  - Bovine, Equine, Porcine, Human, Chicken, Turkey etc.

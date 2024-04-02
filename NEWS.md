# pdtools 0.9.0
  
* Fixes bug in `list_organisms()` where the first organism was being omitted.
* Fixes bug in `download_most_recent_complete()` where the function would be caught in a never ending while loop for organisms with only 1 analysis result.  
* Adds functionality to parse the organisms statistics table: `get_organism_table()` (requires RSelenium and janitor)  


# pdtools 0.8.0

* Fixed issues associated with `.data$` use in tidyselect functions  
* `make_ftp_paths` now warns when genomes are present in the data but not the assembly summary file  
* Updated tests  

# pdtools 0.7.1

* Added a `NEWS.md` file to track changes to the package.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata_cleanup_tools.R
\name{extract_consensus_ag_species}
\alias{extract_consensus_ag_species}
\title{Extract a consensus ag host species from metadata}
\usage{
extract_consensus_ag_species(dat, parallel = FALSE)
}
\arguments{
\item{dat}{an ncbi pathogen detection metadata table}

\item{parallel}{boolean, should furrr be used to parallelize? need to set your future::plan()}
}
\value{
returns a tibble of 2 columns, 1st = target_acc, 2nd = ag_match
}
\description{
Extract a consensus ag host species from metadata
}
\examples{
extract_consensus_ag_species(klebsiella_example_dat)
}

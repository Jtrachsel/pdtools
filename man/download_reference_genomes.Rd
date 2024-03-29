% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome_download_tools.R
\name{download_reference_genomes}
\alias{download_reference_genomes}
\title{Download reference genomes}
\usage{
download_reference_genomes(genome_names, type, data_dir, avail = FALSE)
}
\arguments{
\item{genome_names}{one of LT2, USDA15WA1, FSIS1502916}

\item{type}{fna, gff, etc...}

\item{data_dir}{destination directory}

\item{avail}{TRUE/FALSE, if TRUE, return available genome_name options}
}
\value{
tibble of download results (or vector of genome_names if avail == TRUE)
}
\description{
Download reference genomes
}
\examples{
# download_reference_genome('LT2')
}

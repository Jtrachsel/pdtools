% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_pan_dist}
\alias{example_pan_dist}
\title{A distance matrix describing the jaccard distance between genomes in the
provided example pangenome}
\format{
A dist object
}
\source{
dist(method = 'binary', t(example_pangenome_matrix))
}
\usage{
example_pan_dist
}
\description{
Used for pangenome function examples
}
\keyword{datasets}

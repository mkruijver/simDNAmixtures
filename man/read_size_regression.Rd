% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_size_regression.R
\name{read_size_regression}
\alias{read_size_regression}
\title{Reads a stutter regression file}
\usage{
read_size_regression(filename)
}
\arguments{
\item{filename}{Character. Path to file.}
}
\description{
Reads a stutter regression file
}
\details{
Reads the file from disk and returns a function that provides the fragment length (bp) for a given locus and allele.
}
\examples{
filename <- system.file("extdata","GlobalFiler_SizeRegression.csv",package = "SimMixDNA")
regression <- read_size_regression(filename)
regression("vWA", 12)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiniBenchmark.R
\name{MiniBenchmark}
\alias{MiniBenchmark}
\title{Perform a quick test of the POUMM example likelihood calculation}
\usage{
MiniBenchmark(N = 10000, Ntests = 10)
}
\arguments{
\item{N}{number of tips in the test phylogenetic tree, default 10000}

\item{Ntests}{number of calculations within a call to sys.time (the resulting
times are averages Ntests calls). Default: 10.}
}
\value{
a data.frame containing the response times for the POUMM log-likelihood
calculation.
}
\description{
This function runs a small benchmark to evaluate the SPLITT-package 
installation on a given computer.
}

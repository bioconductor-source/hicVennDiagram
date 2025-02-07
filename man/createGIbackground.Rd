% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gleamTest.R
\name{createGIbackground}
\alias{createGIbackground}
\title{Create background by input GInteractions}
\usage{
createGIbackground(gi, size = 2 * lengths(gi))
}
\arguments{
\item{gi}{A vector of bedpe files or a list of genomic
interaction data (\link[S4Vectors:Pairs-class]{Pairs} or
\link[InteractionSet:GInteractions-class]{GInteractions}).}

\item{size}{The maximal size of the background}
}
\description{
Create background based on the distance distribution of input GInteractions.
}
\examples{

pd <- system.file("extdata", package = "hicVennDiagram", mustWork = TRUE)
fs <- dir(pd, pattern = ".bedpe", full.names = TRUE)[1]
set.seed(123)
# createGIbackground(fs)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitness.R
\name{fitness}
\alias{fitness}
\title{optimization function}
\usage{
fitness(
  k,
  pca.to.tune,
  data.prep,
  IData,
  use = "everything",
  method = "pearson"
)
}
\arguments{
\item{k}{Number of eigenvectors to test}

\item{pca.to.tune}{Pca used for optimization}

\item{data.prep}{object returned by dataPrep function}

\item{IData}{Intensity correlation data}

\item{use}{An optional character string giving a method for computing correlations 
in the presence of missing values. Default is "everything", but when
missing values are present, "pairwise.complete.obs" is required.}
}
\value{
Results for k optimization
}
\description{
Optimization function for k parameter
}
\keyword{internal}

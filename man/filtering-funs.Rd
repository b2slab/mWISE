% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterBased_filter.R, R/clustering_subfuncs.R,
%   R/eps_optimization.R, R/featuresClustering.R, R/k_optimization.R,
%   R/recoveringPeaks.R
\name{filtering-funs}
\alias{filtering-funs}
\alias{clusterBased.filter}
\alias{dataPrep}
\alias{.LaplacianNg}
\alias{eps.optimization}
\alias{featuresClustering}
\alias{k.optimization}
\alias{recoveringPeaks}
\title{Functions to apply cluster-based filtering}
\usage{
clusterBased.filter(
  df,
  Add.Id = NULL,
  Freq = 0.5,
  Info.Add = NULL,
  polarity,
  do.Par = TRUE,
  nClust = 2
)

dataPrep(IData, Rt, Rt.05 = 5, use = "everything", method = "pearson")

eps.optimization(
  pca.to.tune,
  data.prep,
  IData,
  use = "everything",
  k.tuned,
  method = "pearson",
  do.Par,
  nClust
)

featuresClustering(
  Peak.List,
  Intensity.idx,
  use = "everything",
  method = "pearson",
  Rt.05 = 5,
  do.Par = TRUE,
  nClust
)

k.optimization(
  pca.to.tune,
  data.prep,
  IData,
  nrow.List,
  use = "everything",
  method = "pearson",
  do.Par = TRUE,
  nClust
)

recoveringPeaks(Annotated.Tab, MH.Tab)
}
\arguments{
\item{df}{Columns may contain: "Compound", "Add.Id", "Isotope"
"Compound" for the proposed candidates.
"Add.Id" for the adduct or fragment proposed.
"Isotope" to identify the proposed isotopologues.}

\item{Add.Id}{It indicates the adduct(s) or fragment(s) that are 
required to exist.
If NULL, those adducts with an observed frequency equal or 
higher than 0.50 will be used.}

\item{Freq}{Minimum observed frequency to consider an adduct or a fragment 
to apply the filter (Def: 0.5).}

\item{Info.Add}{Data frame with adducts and in source fragments information. 
If NULL, the default mWISE table will be loaded. 
The columns should be:
\itemize{
\item{name} with the name of the adduct or fragment
\item{nmol}{ with the number of molecules (i.e., 2M+H: nmol=2 )}
\item{charge}{ with the charge of the adduct or 
fragment (i.e., M+3H: charge=3)}
\item{massdiff}{ with the mass difference 
(i.e., M+H: massdiff=1.007276)}
\item{quasi}{ with a 1 if the adduct should be considered 
as quasi-molecular and a 0, otherwise (Optional)}
\item{polarity}{ with a character vector indicating the polarity 
of the adduct or fragment. The options are "positive" 
or "negative"}}}

\item{polarity}{Acquisition mode of the study. It can be "positive" 
or "negative".}

\item{do.Par}{TRUE if parallel computing is required. Def: TRUE}

\item{nClust}{Number of cores that may be used if do.Par = TRUE.}

\item{IData}{Data frame containing the intensity for 
each sample in its columns.}

\item{Rt}{Vector containing the retention times.}

\item{Rt.05}{Retention time value to get a similarity of 0.5.}

\item{use}{An optional character string giving a method for computing correlations 
in the presence of missing values. Default is "everything", but when
missing values are present, "pairwise.complete.obs" is required.}

\item{method}{A character string indicating which correlation coefficient
is to be computed. One of "pearson" (default), "kendall", or "spearman".}

\item{pca.to.tune}{PCA to perform the spectral clustering.}

\item{data.prep}{Result returned by \code{dataPrep} function.}

\item{k.tuned}{Optimized number of clusters (k) computed using 
\code{k.optimization} function.}

\item{Peak.List}{Data frame containing the LC-MS features.
Columns should contain:
\itemize{
\item{Peak.Id}{ for a peak identifier}
\item{mz}{ for a mass-to-charge ratio value}
\item{rt}{ for the retention time}
\item{Intensities for each sample}
}}

\item{Intensity.idx}{Numeric vector indicating the column index 
for the intensities}

\item{nrow.List}{Numeric vector indicating the number of peaks.}

\item{Annotated.Tab}{Data frame returned by \code{matchingStage} function.}

\item{MH.Tab}{Data frame returned by \code{clusterBased.filter} function.}
}
\value{
Function \code{clusterBased.filter} returns a data frame of 
filtered candidates.

Function \code{dataPrep} returns a list containing the 
Gaussian similarity 
matrices for the retention time differences
and the intensities correlation.

Function \code{eps.optimization} returns an optimized epsilon 
parameter for dbscan algorithm.

Function \code{featuresClustering} returns the input peak list 
with an additional column named pcgroup
that indicates the clustering.

Function \code{k.optimization} returns the ptimized number of 
clusters (k) using kmeans algorithm.

Function \code{recoveringPeaks} returns a data frame of 
filtered candidates but 
with all peaks recovered.
}
\description{
Function \code{clusterBased.filter} selects only those 
proposed compounds 
for which a quasimolecular adduct or fragment has been also 
proposed in another peak of the same cluster.

Function \code{dataPrep} prepares the intensity and retention 
time data for 
spectral clustering.
Function \code{.LaplacianNg} computes a normalized Laplacian matrix.

Function \code{eps.optimization} optimizes the epsilon 
parameter of the dbscan algorithm.

Function \code{featuresClustering} performs spectral clustering 
to group those features
that come from the same metabolite. It uses \code{dataPrep}, 
\code{.LaplacianNg}, \code{k.optimization}
and \code{eps.optimization} functions. The correlation is computed using the
function \code{cor(use = "pairwise.complete.obs")}.

Function \code{k.optimization} optimizes the number of clusters. 
This value will be used to define the number
of eigenvectors considered in the spectral clustering.

Function \code{recoveringPeaks} recovers the peaks that 
have been removed 
from the first annotated object.
}
\examples{
data("sample.keggDB")
Cpd.Add <- CpdaddPreparation(KeggDB = sample.keggDB, do.Par = FALSE)
data(sample.dataset)
Peak.List <- sample.dataset$Positive$Input
Annotated.List <- matchingStage(Peak.List = Peak.List, Cpd.Add = Cpd.Add,
                                polarity = "positive", do.Par = FALSE)
Intensity.idx <- seq(27,38)
clustered <- featuresClustering(Peak.List = Peak.List, 
                                Intensity.idx = Intensity.idx, 
                                do.Par = FALSE)
Annotated.Tab <- Annotated.List$Peak.Cpd
Annotated.Tab <- merge(Annotated.Tab,
                       clustered$Peak.List[,c("Peak.Id", "pcgroup")],
                       by = "Peak.Id")
                       
MH.Tab <- clusterBased.filter(df = Annotated.Tab, 
                              polarity = "positive")
                              
data(sample.dataset)
Peak.List <- sample.dataset$Positive$Input
Intensity.idx <- seq(27,38)
clustered <- featuresClustering(Peak.List = Peak.List, 
                                Intensity.idx = Intensity.idx, 
                                do.Par = FALSE)
data("sample.keggDB")
Cpd.Add <- CpdaddPreparation(KeggDB = sample.keggDB, 
do.Par = FALSE)
data(sample.dataset)
Peak.List <- sample.dataset$Positive$Input
Annotated.List <- matchingStage(Peak.List = Peak.List, Cpd.Add = Cpd.Add,
                                polarity = "positive", do.Par = FALSE)
Intensity.idx <- seq(27,38)
clustered <- featuresClustering(Peak.List = Peak.List, 
                                Intensity.idx = Intensity.idx, 
                                do.Par = FALSE)
Annotated.Tab <- Annotated.List$Peak.Cpd
Annotated.Tab <- merge(Annotated.Tab,
                       clustered$Peak.List[,c("Peak.Id", "pcgroup")],
                       by = "Peak.Id")
                       
MH.Tab <- clusterBased.filter(df = Annotated.Tab, 
                              polarity = "positive")
                              
recoveredPeaks <- recoveringPeaks(Annotated.Tab = Annotated.Tab,
                                  MH.Tab = MH.Tab)
}

#' @name filtering-funs
#' @aliases featuresClustering
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{featuresClustering} performs spectral clustering 
#' to group those features
#' that come from the same metabolite. It uses \code{dataPrep}, 
#' \code{.LaplacianNg}, \code{k.optimization}
#' and \code{eps.optimization} functions. The correlation is computed using the
#' function \code{cor(use = "pairwise.complete.obs")}.
#' @param Peak.List
#' Data frame containing the LC-MS features.
#' Columns should contain:
#' \itemize{
#' \item{Peak.Id}{ for a peak identifier}
#' \item{mz}{ for a mass-to-charge ratio value}
#' \item{rt}{ for the retention time}
#' \item{Intensities for each sample}
#' }
#' @param Intensity.idx
#' Numeric vector indicating the column index 
#' for the intensities
#' @param Rt.05
#' Retention time value to get a similarity of 0.5.
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of cores that may be used if do.Par = TRUE.
#' @return
#' Function \code{featuresClustering} returns the input peak list 
#' with an additional column named pcgroup
#' that indicates the clustering.
#' @examples 
#' data(sample.dataset)
#' Peak.List <- sample.dataset$Positive$Input
#' Intensity.idx <- seq(27,38)
#' clustered <- featuresClustering(Peak.List = Peak.List, 
#'                                 Intensity.idx = Intensity.idx, 
#'                                 do.Par = FALSE)
#' @export
#' @importFrom dbscan dbscan
#' @importFrom stats prcomp

featuresClustering <- function(Peak.List, Intensity.idx,
                               Rt.05 = 5, do.Par = TRUE, nClust) {

  cat("Preparing the data for clustering...")
  IData <- Peak.List[,Intensity.idx]
  Rt <- Peak.List$rt
  data.prep <- dataPrep(IData = IData, Rt = Rt, Rt.05 = Rt.05)
  comb.mat.clustering <- sqrt(data.prep$Rt.sim*data.prep$I.sim)
  data.prep$Rt.sim <- data.prep$Rt.sim-diag(1,dim(data.prep$Rt.sim)[1],
                                            dim(data.prep$Rt.sim)[2])
  data.prep$I.sim <- data.prep$I.sim-diag(1,dim(data.prep$I.sim)[1],
                                          dim(data.prep$I.sim)[2])
  data.prep$Rt.sim[data.prep$Rt.sim == 0] <- 1e-16
  data.prep$I.sim[data.prep$I.sim == 0] <- 1e-16
  comb.mat <- sqrt(data.prep$Rt.sim*data.prep$I.sim)
  Lapl.mat <- .LaplacianNg(mat = comb.mat)
  pca.tune <- prcomp(Lapl.mat, center = TRUE, scale. = TRUE)
  rownames(pca.tune$x) <- rownames(comb.mat)
  cat("DONE!","\n")
  # Parameters optimization
  cat("Computing optimized parameters for spectral clustering...")
  k.tuned <- k.optimization(pca.to.tune = pca.tune, data.prep = data.prep, 
                            IData = IData, nrow.List = nrow(Peak.List), 
                            do.Par = do.Par, nClust = nClust)
  eps.tuned <- eps.optimization(pca.to.tune = pca.tune, 
                                data.prep = data.prep, 
                                IData = IData, k.tuned = k.tuned, 
                                do.Par = do.Par, nClust = nClust)
  cat("DONE!","\n")
  cat("Clustering peaks...")
  Lapl.mat <- .LaplacianNg(mat = comb.mat.clustering)
  pca.obj <- prcomp(Lapl.mat, center = TRUE, scale. = TRUE)
  rownames(pca.obj$x) <- rownames(comb.mat.clustering)
  # Clustering with tuned parameters
  dbscan.clustering <- dbscan::dbscan(x = pca.obj$x[,seq_len(k.tuned)], 
                                      eps = eps.tuned, minPts = 1)
  Peak.List$pcgroup <- dbscan.clustering$cluster
  cat("DONE!","\n")
  return(list(Peak.List = Peak.List, 
              k.tuned = k.tuned, 
              eps.tuned = eps.tuned))
}

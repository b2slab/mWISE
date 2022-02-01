#' @name filtering-funs
#' @aliases k.optimization
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{k.optimization} optimizes the number of clusters. 
#' This value will be used to define the number
#' of eigenvectors considered in the spectral clustering.
#' @param pca.to.tune
#' PCA to perform the spectral clustering.
#' @param data.prep
#' Result returned by \code{dataPrep} function.
#' @param IData
#' Data frame containing the intensity for 
#' each sample in its columns.
#' @param nrow.List
#' Numeric vector indicating the number of peaks.
#' @param use
#' An optional character string giving a method for computing correlations 
#' in the presence of missing values. Default is "everything", but when
#' missing values are present, "pairwise.complete.obs" is required. 
#' @param method
#' A character string indicating which correlation coefficient
#' is to be computed. One of "pearson" (default), "kendall", or "spearman".
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of cores that may be used if do.Par = TRUE.
#' @return
#' Function \code{k.optimization} returns the ptimized number of 
#' clusters (k) using kmeans algorithm.
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply

k.optimization <- function(pca.to.tune, data.prep, IData, nrow.List, 
                           use = "everything", method = "pearson",
                           do.Par=TRUE, nClust) {
  if (do.Par){
    doParallel::registerDoParallel(nClust)
  }

  ks <- seq(2,(nrow.List/4), 20)
  res.opt <- plyr::ldply(ks, function(k){
    res.opt <- fitness(k = k, pca.to.tune = pca.to.tune, use = use,
                       data.prep = data.prep, IData = IData,
                       method = method)
    return(data.frame(k = k, t(res.opt)))
  },.parallel = do.Par)
  res.opt$X1 <- res.opt$X1/(max(res.opt$X1))
  res.opt$X2 <- res.opt$X2/(max(res.opt$X2))
  res.opt$X3 <- res.opt$X3/(max(res.opt$X3))
  res.opt$s <-res.opt$X3-((res.opt$X1+res.opt$X2)/2)
  k.tuned <- res.opt$k[which.min(res.opt$s)]
  return(k.tuned)

}


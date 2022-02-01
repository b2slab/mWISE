#' @name filtering-funs
#' @aliases eps.optimization
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{eps.optimization} optimizes the epsilon 
#' parameter of the dbscan algorithm.
#' @param pca.to.tune
#' PCA to perform the spectral clustering.
#' @param data.prep
#' Result returned by \code{dataPrep} function.
#' @param IData
#' Data frame containing the intensity for each 
#' sample in its columns.
#' @param k.tuned
#' Optimized number of clusters (k) computed using 
#' \code{k.optimization} function.
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
#' Function \code{eps.optimization} returns an optimized epsilon 
#' parameter for dbscan algorithm.
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply
#' @importFrom stats dist
#' 
eps.optimization <- function(pca.to.tune, data.prep, IData, 
                             use = "everything", k.tuned, 
                             method = "pearson", do.Par, nClust) {
  if (do.Par){
    doParallel::registerDoParallel(nClust)
  }
  dist.cl.rand <- plyr::ldply(seq_len(1000), function(i){
    l <- sample(2:10,1)
    M <- mean(dist(pca.to.tune$x[sample(seq_len(nrow(pca.to.tune$x)), 
                                        size = l),
                                 seq_len(k.tuned)]))
    data.frame(clust = i, Size = l, Mean = M)
  })
  min.value <- mean(dist.cl.rand$Mean)/2

  eps.opt <- seq(0.1,round(min.value, digits = 1),by = 0.1)
  res.opt <- plyr::ldply(eps.opt, function(e){
    res.opt <- fitness.eps(eps = e, pca.to.tune = pca.to.tune, 
                           k.tuned = k.tuned, IData = IData, 
                           data.prep = data.prep, use = use,
                           method = method)
    return(data.frame(k = e, t(res.opt)))
  },.parallel = do.Par)

  res.opt$X1 <- res.opt$X1/(max(res.opt$X1))
  res.opt$X2 <- res.opt$X2/(max(res.opt$X2))
  res.opt$X3 <- res.opt$X3/(max(res.opt$X3))

  res.opt$s <-res.opt$X3-(res.opt$X1+res.opt$X2)
  eps.tuned <- res.opt$k[which.min(res.opt$s)]
  return(eps.tuned)
}

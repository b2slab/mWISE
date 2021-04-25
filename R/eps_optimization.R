#' @name filtering-funs
#' @aliases k.optimization
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{eps.optimization} optimizes the epsilon parameter of the dbscan algorithm.
#' @example
#' eps <- eps.optimization(pca.to.tune, data.prep, IData, k.tuned, nrow.List, do.Par, nClust)
#' @param pca.to.tune
#' PCA to perform the spectral clustering.
#' @param data.prep
#' Result returned by \code{dataPrep} function.
#' @param IData
#' Data frame containing the intensity for each sample in its columns.
#' @param k.tuned
#' Optimized number of clusters (k) computed using \code{k.optimization} function.
#' @param nrow.List
#' Numeric vector indicating the number of peaks.
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of cores that may be used if do.Par = TRUE.
#' @return
#' Function \code{eps.optimization} returns an optimized epsilon parameter for dbscan algorithm.
#' @export

eps.optimization <- function(pca.to.tune, data.prep, IData, k.tuned, do.Par, nClust) {
  doParallel::registerDoParallel(nClust)
  fitness.eps <- function(eps, pca.to.tune, k.tuned){
    dbscan.clustering <- dbscan::dbscan(x = pca.to.tune$x[,1:k.tuned], eps = eps)

    # calculate the I similarity mean of each cluster
    Mean.I <- plyr::llply(unique(dbscan.clustering$cluster), function(c){
      mean(data.prep$I.sim[which(dbscan.clustering$cluster %in% c),which(dbscan.clustering$cluster %in% c)])
    }) #%>% unlist() %>% sum()
    Mean.I <- sum(unlist(Mean.I))

    # calculate the RT similarity mean of each cluster
    Mean.RT <- plyr::llply(unique(dbscan.clustering$cluster), function(c){
      mean(data.prep$Rt.sim[which(dbscan.clustering$cluster %in% c),which(dbscan.clustering$cluster %in% c)])
    }) #%>% unlist() %>% sum()
    Mean.RT <- sum(unlist(Mean.RT))

    n.Corr <- plyr::ldply(unique(dbscan.clustering$cluster), function(c){
      IData[which.max(rowMeans(IData[which(dbscan.clustering$cluster %in% c),])),]
    })
    #cor.mat <- n.Corr %>% t() %>% cor()
    cor.mat <- cor(t(n.Corr))
    cor.mat[upper.tri(x = cor.mat, diag = T)] <- 0
    n <- sum(cor.mat>0.6)

    return(c(Mean.I, Mean.RT, n))

  }

  dist.cl.rand <- plyr::ldply(1:1000, function(i){
    l <- sample(2:10,1)
    M <- mean(dist(pca.to.tune$x[sample(1:nrow(pca.to.tune$x), size = l),1:k.tuned]))
    data.frame(clust = i, Size = l, Mean = M)
  })
  min.value <- mean(dist.cl.rand$Mean)/2

  eps.opt <- seq(0.1,round(min.value, digits = 1),by = 0.1)
  res.opt <- plyr::ldply(eps.opt, function(e){
    res.opt <- fitness.eps(eps = e, pca.to.tune = pca.to.tune, k.tuned = k.tuned)
    return(data.frame(k = e, t(res.opt)))
  },.parallel = do.Par)

  res.opt$X1 <- res.opt$X1/(max(res.opt$X1))
  res.opt$X2 <- res.opt$X2/(max(res.opt$X2))
  res.opt$X3 <- res.opt$X3/(max(res.opt$X3))

  res.opt$s <-res.opt$X3-(res.opt$X1+res.opt$X2)
  eps.tuned <- res.opt$k[which.min(res.opt$s)]
  return(eps.tuned)
}

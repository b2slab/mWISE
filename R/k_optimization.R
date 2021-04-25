#' @name filtering-funs
#' @aliases k.optimization
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{k.optimization} optimizes the number of clusters. This value will be used to define the number
#' of eigenvectors considered in the spectral clustering.
#' @example
#' k <- k.optimization(pca.to.tune, data.prep, IData, nrow.List, do.Par, nClust)
#' @param pca.to.tune
#' PCA to perform the spectral clustering.
#' @param data.prep
#' Result returned by \code{dataPrep} function.
#' @param IData
#' Data frame containing the intensity for each sample in its columns.
#' @param nrow.List
#' Numeric vector indicating the number of peaks.
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of cores that may be used if do.Par = TRUE.
#' @return
#' Function \code{k.optimization} returns the ptimized number of clusters (k) using kmeans algorithm.
#' @export

k.optimization <- function(pca.to.tune, data.prep, IData, nrow.List, do.Par, nClust) {
  doParallel::registerDoParallel(nClust)
  fitness <- function(k, pca.to.tune){
    kmeans.clustering <- kmeans(x = pca.to.tune$x[,1:k], centers = k)

    # calculate the I similarity mean of each cluster
    Mean.I <- plyr::llply(unique(kmeans.clustering$cluster), function(c){
      mean(data.prep$I.sim[which(kmeans.clustering$cluster %in% c),which(kmeans.clustering$cluster %in% c)])
    })# %>% unlist() %>% sum()
    Mean.I <- sum(unlist(Mean.I))

    # calculate the RT similarity mean of each cluster
    Mean.RT <- plyr::llply(unique(kmeans.clustering$cluster), function(c){
      mean(data.prep$Rt.sim[which(kmeans.clustering$cluster %in% c),which(kmeans.clustering$cluster %in% c)])
    }) #%>% unlist() %>% sum()
    Mean.RT <- sum(unlist(Mean.RT))

    # Number of putative compound units correlated
    n.Corr <- plyr::ldply(unique(kmeans.clustering$cluster), function(c){
      colMeans(IData[which(kmeans.clustering$cluster %in% c),])
    })
    #cor.mat <- n.Corr %>% t() %>% cor()
    cor.mat <- cor(t(n.Corr))
    cor.mat[upper.tri(x = cor.mat, diag = T)] <- 0
    n <- sum(cor.mat>0.6)

    return(c(Mean.I, Mean.RT, n))

  }

  ks <- seq(2,(nrow.List/4), 20)
  res.opt <- plyr::ldply(ks, function(k){
    res.opt <- fitness(k = k, pca.to.tune = pca.to.tune)
    return(data.frame(k = k, t(res.opt)))
  },.parallel = do.Par)
  res.opt$X1 <- res.opt$X1/(max(res.opt$X1))
  res.opt$X2 <- res.opt$X2/(max(res.opt$X2))
  res.opt$X3 <- res.opt$X3/(max(res.opt$X3))
  res.opt$s <-res.opt$X3-((res.opt$X1+res.opt$X2)/2)
  k.tuned <- res.opt$k[which.min(res.opt$s)]
  return(k.tuned)

}


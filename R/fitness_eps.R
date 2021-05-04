#' @title optimization function
#' @description  
#' Optimization function for epsilon parameter
#' @param eps
#' Epsilon parameter test
#' @param pca.to.tune
#' Pca used for optimization
#' @param k.tuned
#' Optimized k value
#' @param data.prep
#' object returned by dataPrep function
#' @param IData
#' Intensity correlation data
#' @return 
#' results for epsilon optimization
#' @keywords internal
#' @importFrom dbscan dbscan
#' @importFrom plyr ldply llply
#' @importFrom stats cor
fitness.eps <- function(eps, pca.to.tune, k.tuned, data.prep, IData){
  dbscan.clustering <- dbscan::dbscan(x = pca.to.tune$x[,seq_len(k.tuned)], 
                                      eps = eps)
  
  # calculate the I similarity mean of each cluster
  Mean.I <- plyr::llply(unique(dbscan.clustering$cluster), function(c){
    mean(data.prep$I.sim[which(dbscan.clustering$cluster %in% c),
                         which(dbscan.clustering$cluster %in% c)])
  }) #%>% unlist() %>% sum()
  Mean.I <- sum(unlist(Mean.I))
  
  # calculate the RT similarity mean of each cluster
  Mean.RT <- plyr::llply(unique(dbscan.clustering$cluster), function(c){
    mean(data.prep$Rt.sim[which(dbscan.clustering$cluster %in% c),
                          which(dbscan.clustering$cluster %in% c)])
  }) #%>% unlist() %>% sum()
  Mean.RT <- sum(unlist(Mean.RT))
  
  n.Corr <- plyr::ldply(unique(dbscan.clustering$cluster), function(c){
    IData[which.max(
      rowMeans(IData[which(dbscan.clustering$cluster %in% c),])),]
  })
  #cor.mat <- n.Corr %>% t() %>% cor()
  cor.mat <- cor(t(n.Corr))
  cor.mat[upper.tri(x = cor.mat, diag = TRUE)] <- 0
  n <- sum(cor.mat>0.6)
  
  return(c(Mean.I, Mean.RT, n))
  
}
#' @title optimization function
#' @description  
#' Optimization function for k parameter
#' @param k
#' Number of eigenvectors to test
#' @param pca.to.tune
#' Pca used for optimization
#' @param data.prep
#' object returned by dataPrep function
#' @param IData
#' Intensity correlation data
#' @param use
#' An optional character string giving a method for computing correlations 
#' in the presence of missing values. Default is "everything", but when
#' missing values are present, "pairwise.complete.obs" is required. 
#' @return 
#' Results for k optimization
#' @keywords internal
#' @importFrom stats kmeans cor
#' @importFrom plyr ldply llply

fitness <- function(k, pca.to.tune, data.prep, IData, use = "everything",
                    method = "pearson"){
  kmeans.clustering <- kmeans(x = pca.to.tune$x[,seq_len(k)], centers = k)
  
  # calculate the I similarity mean of each cluster
  Mean.I <- plyr::llply(unique(kmeans.clustering$cluster), function(c){
    mean(data.prep$I.sim[which(kmeans.clustering$cluster %in% c),
                         which(kmeans.clustering$cluster %in% c)])
  })# %>% unlist() %>% sum()
  Mean.I <- sum(unlist(Mean.I))
  
  # calculate the RT similarity mean of each cluster
  Mean.RT <- plyr::llply(unique(kmeans.clustering$cluster), function(c){
    mean(data.prep$Rt.sim[which(kmeans.clustering$cluster %in% c),
                          which(kmeans.clustering$cluster %in% c)])
  }) #%>% unlist() %>% sum()
  Mean.RT <- sum(unlist(Mean.RT))
  
  # Number of putative compound units correlated
  n.Corr <- plyr::ldply(unique(kmeans.clustering$cluster), function(c){
    colMeans(IData[which(kmeans.clustering$cluster %in% c),], na.rm = TRUE)
  })
  #cor.mat <- n.Corr %>% t() %>% cor()
  cor.mat <- cor(t(n.Corr), use = use, method = method)
  cor.mat[upper.tri(x = cor.mat, diag = TRUE)] <- 0
  n <- sum(cor.mat>0.6)
  return(c(Mean.I, Mean.RT, n))
}

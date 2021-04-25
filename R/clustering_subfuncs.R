#' @name filtering-funs
#' @aliases dataPrep
#' @aliases .LaplacianNg
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{dataPrep} prepares the intensity and retention time data for spectral clustering.
#' Function \code{.LaplacianNg} computes a normalized Laplacian matrix.
#' @param IData
#' Data frame containing the intensity for each sample in its columns.
#' @param Rt
#' Vector containing the retention times.
#' @param Rt.05
#' Retention time value to get a similarity of 0.5.
#' @param mat
#' Matrix to compute the normalized Laplacian matrix.
#' @return
#' Function \code{dataPrep} returns a list containing the Gaussian similarity matrices for the retention time differences
#' and the intensities correlation.
#' @export

# data preparation
dataPrep <- function(IData, Rt, Rt.05 = 5) {
  I.sim <- cor(t(IData))
  I.dist <- (1-I.sim)/2
  I.05 <- 0.4
  I.05.dist <- (1-I.05)/2
  sigma.I <- sqrt(-1*((I.05.dist^2)/(2*log(0.5))))
  I <- exp(-(((I.dist)^2)/(2*(sigma.I^2))))
  # 0's are dangerous afterwards, they should be replaced by sth safer
  I[I == 0] <- 1e-16
  # Rt matrix
  RtData <- abs(outer(Rt, Rt, '-'))
  rownames(RtData) <- rownames(I)
  colnames(RtData) <- colnames(I)
  sigma.rt <- sqrt(-1*((Rt.05^2)/(2*log(0.5))))
  Rt <- exp(-(((RtData)^2)/(2*(sigma.rt^2))))
  # 0's are dangerous afterwards, they should be replaced by sth safer
  Rt[Rt == 0] <- 1e-16
  return(list(I.sim = I, Rt.sim = Rt))
}


# Laplacian kernel
.LaplacianNg <- function(mat)
{
  D <- rowSums(mat)
  sqriD <- diag(1/sqrt(D))
  return(sqriD %*% mat %*% sqriD)
}

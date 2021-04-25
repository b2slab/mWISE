#' @name matching-funs
#' @aliases matchingStage
#' @title Functions to match a peak-intensity table to KEGG database
#' @description
#' Function \code{mz.ppm.converter} transforms the user-defined tolerance from ppm to mz uncertainty or viceversa.
#' @param tol
#' Tolerance in ppm. Def: 10
#' @param mzrange
#' mz uncertainty. Def: 0.001
#' @param mz
#' mass-to-charge ratio for which the conversion is performed.
#' @param conversion
#' Direction of the conversion. Character vector containing some of: "ppm.to.mz", "mz.to.ppm"
#' @return
#' Function \code{mz.ppm.converter} returns a vector of values containing uncertainty (if conversion="ppm.to.mz")
#' or tolerance in ppm (if conversion="mz.to.ppm")
#' @export

mz.ppm.converter <- function(tol = 10, mzrange = 0.001, mz,
                             conversion = c("ppm.to.mz", "mz.to.ppm")){
  if (conversion == "ppm.to.mz"){
    uncertainty <- mz*tol*10^(-6)
    return(uncertainty)
  }
  if (conversion == "mz.to.ppm"){
    tol <- mzrange/(mz*10^(-6))
    return(tol)
  }
}

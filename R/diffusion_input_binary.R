#' @name diffusion-funs
#' @aliases diffusion.input.Binary
#' @title Functions to apply diffusion in graphs
#' @description
#' Function \code{diffusion.input.Binary} computes the binary diffusion input score.
#' @example
#' Diffusion.Input<-diffusion.input.Binary(df = Peak.Cpd,do.Par = T,  nClust = detectCores()-1)
#' @param df
#' Data frame containing the potential candidates. It is recommended to use the
#' data frame resulted from mWISE clustered-based filtering. Columns may contain "Peak.Id" for a peak identifier and
#' "Compound" for a KEGG ID.
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{diffusion.input.Binary} returns a data frame containing the binary diffusion input.
#' @export

diffusion.input.Binary<-function(df, do.Par = TRUE, nClust = 2){
  cat("Computing diffusion input...")
  cat("DONE!","\n")
  return(data.frame(Compound = unique(df$Compound),
                    Diffusion.Input=1))
}

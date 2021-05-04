#' @name diffusion-funs
#' @aliases diffusion.input.probability
#' @title Functions to apply diffusion in graphs
#' @description
#' Function \code{diffusion.input.probability} computes the 
#' probability diffusion input score.
#' @param df
#' Data frame containing the potential candidates. 
#' It is recommended to use the
#' data frame resulted from mWISE clustered-based filtering. 
#' Columns may contain "Peak.Id" for a peak identifier and
#' "Compound" for a KEGG ID.
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{diffusion.input.probability} returns a data 
#' frame containing the probability diffusion input.
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply
diffusion.input.probability <- function(df, do.Par = TRUE, nClust = 2) {
  if (do.Par){
    doParallel::registerDoParallel(nClust)
  }
  cat("Computing diffusion input...")
  Prob.Tabs <- plyr::ldply(as.character(unique(df$Compound)), function(cpd){
    dx <- df[df$Compound %in% cpd,]
    p.peaks <- plyr::llply(dx$Peak.Id, function(id) {
      dy <- df[df$Peak.Id %in% id,]
      Ppk <- 1/nrow(dy)
      return(1-Ppk)
    })# %>% unlist()
    p.peaks <- unlist(p.peaks)
    Pci <- 1-prod(p.peaks)
    data.frame(Compound = cpd, Diffusion.Input = Pci)
  },.parallel = do.Par)
  cat("DONE!","\n")
  return(Prob.Tabs)
}

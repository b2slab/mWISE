#' @name CpdaddPreparation
#' @title Function to compute a compound to adduct matrix
#' @description
#' Function \code{CpdaddPreparation} computes a matrix that allows the fast matching of a peak list 
#' to KEGG database from a data frame containing adducts and fragments information.
#' @example
#' CpdAdd <- CpdaddPreparation(Info.Add, do.Par = TRUE, nClust)
#' @param Info.Add
#' Data frame with adducts and in source fragments information. 
#' The columns should be:
#' \itemize{
#' \item{name} with the name of the adduct or fragment
#' \item{nmol}{ with the number of molecules (i.e., 2M+H: nmol=2 )}
#' \item{charge}{ with the charge of the adduct or fragment (i.e., M+3H: charge=3)}
#' \item{massdiff}{ with the mass difference (i.e., M+H: massdiff=1.007276)}
#' \item{quasi}{ with a 1 if the adduct should be considered as quasi-molecular and a 0, otherwise (Optional)}
#' \item{polarity}{ with a character vector indicating the polarity of the adduct or fragment. The options are "positive" 
#' or "negative"}
#' }
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE.
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{CpdaddPreparation} returns a data frame with the mass-to-charge ratio of all the KEGG compounds and
#' adducts or fragments provided in Info.Add.
#' @export

CpdaddPreparation <- function(Info.Add, do.Par = TRUE, nClust) {
  doParallel::registerDoParallel(nClust)
  Cpd.Add.cmps <- mWISE::Cpd.Add[,c("Compound", "exact_mass")]
  Cpd.Add.cmps <- Cpd.Add.cmps[!duplicated(Cpd.Add.cmps$Compound),]
  New.Cpd.Add <- plyr::ldply(1:nrow(Cpd.Add.cmps), function(i){
    dx <- Cpd.Add.cmps[i,]
    adds.mz <- ((dx$exact_mass*Info.Add$nmol)+Info.Add$massdiff)/abs(Info.Add$charge)
    cpd.add <- data.frame(Compound = dx$Compound, exact_mass = dx$exact_mass, 
                          Add.name = Info.Add$name, Polarity = Info.Add$polarity, 
                          mz.Add = adds.mz)
    return(cpd.add)
  },.parallel = do.Par)
}
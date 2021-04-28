#' @name adducts-funs
#' @aliases printAdducts
#' @aliases printQM_Adducts
#' @title Functions to assist adducts and fragments selection
#' @description
#' Function \code{printAdducts} returns a list of all the available adducts and fragments in mWISE for a given polarity.
#' Function \code{printQM_Adducts} returns a list of the default quasi-molecular adducts or fragments of mWISE for a given polarity.
#' These functions can be used to ease the selection of adducts or fragments to be considered during the matching stage.
#' @param pol
#' Acquisition mode. It can be "positive", "negative" or both. 
#' @export
printAdducts <- function(pol = c("positive", "negative")){
  Info.Add <- mWISE::Info.Add
  Adds <- as.character(Info.Add$name[Info.Add$polarity %in% pol])
  return(Adds)
}
printQM_Adducts <- function(pol = c("positive", "negative")){
  Info.Add <- mWISE::Info.Add
  QM.InfoAdd <- Info.Add[Info.Add$quasi %in% 1,]
  QM.Adds <- as.character(QM.InfoAdd$name[QM.InfoAdd$polarity %in% pol])
  return(QM.Adds)
}

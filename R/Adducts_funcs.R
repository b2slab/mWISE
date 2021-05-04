#' @name adducts-funs
#' @aliases printAdducts
#' @title Functions to assist adducts and fragments selection
#' @description
#' Function \code{printAdducts} returns a list of all the available adducts 
#' and fragments in mWISE for a given polarity.
#' @param pol
#' Acquisition mode. It can be positive, negative or both. 
#' @return
#' Function \code{printAdducts} returns the complete list of adducts 
#' available in mWISE for a given polarity.
#' @examples
#' Adducts.List <- printAdducts(pol="negative")
#' @export
#' @importFrom utils data

printAdducts <- function(pol = c("positive", "negative")){
    Info.Add <- NULL
    data("Info.Add", envir=environment())
    Adds <- as.character(Info.Add$name[Info.Add$polarity %in% pol])
    return(Adds)
}
#' @name adducts-funs
#' @aliases printQM_Adducts
#' @title Functions to assist adducts and fragments selection
#' @description
#' Function \code{printQM_Adducts} returns a list of the default 
#' quasi-molecular adducts or fragments of mWISE for a given polarity.
#' These functions can be used to ease the selection of adducts 
#' or fragments to be considered during the matching stage.
#' @param pol
#' Acquisition mode. It can be positive, negative or both. 
#' @return
#' Function \code{printQM_Adducts} returns the complete list of 
#' quasi-molecular adducts available in mWISE for a given polarity.
#' @examples
#' Adducts.List <- printQM_Adducts(pol="negative")
#' @export
#' @importFrom utils data
printQM_Adducts <- function(pol = c("positive", "negative")){
    Info.Add <- NULL
    data("Info.Add", envir=environment())
    QM.InfoAdd <- Info.Add[Info.Add$quasi %in% 1,]
    QM.Adds <- as.character(QM.InfoAdd$name[QM.InfoAdd$polarity %in% pol])
    return(QM.Adds)
}

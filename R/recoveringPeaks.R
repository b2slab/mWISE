#' @name filtering-funs
#' @aliases recoveringPeaks
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{recoveringPeaks} recovers the peaks that have been removed 
#' from the first annotated object.
#' @param Annotated.Tab
#' Data frame returned by \code{matchingStage} function.
#' @param MH.Tab
#' Data frame returned by \code{clusterBased.filter} function.
#' @return
#' Function \code{recoveringPeaks} returns a data frame of filtered candidates but 
#' with all peaks recovered.
#' @export

recoveringPeaks <- function(Annotated.Tab, MH.Tab) {
  peak.ids <- unique(Annotated.Tab$Peak.Id)[!(unique(Annotated.Tab$Peak.Id) %in% unique(MH.Tab$Peak.Id))]
  recoveredTab <- Annotated.Tab[Annotated.Tab$Peak.Id %in% peak.ids,]
  rbind(MH.Tab, recoveredTab)
}
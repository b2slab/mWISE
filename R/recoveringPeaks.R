#' @name filtering-funs
#' @aliases recoveringPeaks
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{recoveringPeaks} recovers the peaks that 
#' have been removed 
#' from the first annotated object.
#' @param Annotated.Tab
#' Data frame returned by \code{matchingStage} function.
#' @param MH.Tab
#' Data frame returned by \code{clusterBased.filter} function.
#' @return
#' Function \code{recoveringPeaks} returns a data frame of 
#' filtered candidates but 
#' with all peaks recovered.
#' @examples 
#' data("sample.keggDB")
#' Cpd.Add <- CpdaddPreparation(KeggDB = sample.keggDB, 
#' do.Par = FALSE)
#' data(sample.dataset)
#' Peak.List <- sample.dataset$Positive$Input
#' Annotated.List <- matchingStage(Peak.List = Peak.List, Cpd.Add = Cpd.Add,
#'                                 polarity = "positive", do.Par = FALSE)
#' Intensity.idx <- seq(27,38)
#' clustered <- featuresClustering(Peak.List = Peak.List, 
#'                                 Intensity.idx = Intensity.idx, 
#'                                 do.Par = FALSE)
#' Annotated.Tab <- Annotated.List$Peak.Cpd
#' Annotated.Tab <- merge(Annotated.Tab,
#'                        clustered$Peak.List[,c("Peak.Id", "pcgroup")],
#'                        by = "Peak.Id")
#'                        
#' MH.Tab <- clusterBased.filter(df = Annotated.Tab, 
#'                               polarity = "positive")
#'                               
#' recoveredPeaks <- recoveringPeaks(Annotated.Tab = Annotated.Tab,
#'                                   MH.Tab = MH.Tab)
#' @export

recoveringPeaks <- function(Annotated.Tab, MH.Tab) {
  peak.ids <- unique(
    Annotated.Tab$Peak.Id)[!(
      unique(Annotated.Tab$Peak.Id) %in% unique(MH.Tab$Peak.Id))]
  recoveredTab <- Annotated.Tab[Annotated.Tab$Peak.Id %in% peak.ids,]
  rbind(MH.Tab, recoveredTab)
}

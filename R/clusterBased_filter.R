#' @name filtering-funs
#' @aliases clusterBased.filter
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{clusterBased.filter} selects only those proposed compounds for which a quasimolecular
#' adduct or fragment has been also proposed in another peak of the same cluster.
#' @example
#' Filtered.Table <- clusterBased.filter(df = Peak.Cpd, Add.Id = NULL, Freq = 0.10, polarity = "positive")
#' @param df
#' Columns may contain: "Compound", "Add.Id", "Isotope"
#'  "Compound" for the proposed candidates.
#'  "Add.Id" for the adduct or fragment proposed.
#'  "Isotope" to identify the proposed isotopologues.
#' @param Freq
#' Minimum observed frequency to consider an adduct or a fragment to apply the filter.
#' @param Add.Id
#' It indicates the adduct(s) or fragment(s) that are required to exist.
#' If NULL, those adducts with an observed frequency equal or higher than 0.10 will be used.
#' @param polarity
#' Acquisition mode of the study. It can be "positive" or "negative".
#' @return
#' Function \code{clusterBased.filter} returns a data frame of filtered candidates.
#' @export

clusterBased.filter <- function(df, Add.Id = NULL, Freq = 0.50, polarity) {
  if (is.null(Add.Id)){
    Info.Add <- mWISE::Info.Add
    Add.Id.freq <- as.character(Info.Add$name[(Info.Add$Freq>=Freq)&(Info.Add$polarity %in% polarity)])
    Add.Id.QM <- as.character(Info.Add$name[(Info.Add$quasi==1)&(Info.Add$polarity %in% polarity)])
    Add.Id <- unique(c(Add.Id.freq, Add.Id.QM))
  }
  cat(paste("Applying the cluster-based filter using ", paste(Add.Id, collapse = ", "), "...", sep = ""))
  MH.df <- plyr::ldply(unique(df$pcgroup), function(id.group){
    group.df <- df[df$pcgroup %in% id.group,]
    valid.cmps <- as.character(group.df$Compound[(group.df$Add.name %in% Add.Id)])
    return(group.df[group.df$Compound %in% valid.cmps,])
  },.parallel = TRUE)
  cat("DONE!","\n")
  return(MH.df)
}

recoveringPeaks <- function(Annotated.Tab, MH.Tab) {
  peak.ids <- unique(Annotated.Tab$Peak.Id)[!(unique(Annotated.Tab$Peak.Id) %in% unique(MH.Tab$Peak.Id))]
  recoveredTab <- Annotated.Tab[Annotated.Tab$Peak.Id %in% peak.ids,]
  rbind(MH.Tab, recoveredTab)
}

# recoveringPeaks <- function(df, Add.Id = NULL, Freq = 0.10, polarity) {
#   if (is.null(Add.Id)){
#     Info.Add <- mWISE::Info.Add
#     Add.Id.freq <- as.character(Info.Add$name[(Info.Add$Freq>=Freq)&(Info.Add$polarity %in% polarity)])
#     Add.Id.QM <- as.character(Info.Add$name[(Info.Add$quasi==1)&(Info.Add$polarity %in% polarity)])
#     Add.Id <- unique(c(Add.Id.freq, Add.Id.QM))
#   }
#   MH.df <- plyr::ldply(unique(df$pcgroup), function(id.group){
#     group.df <- df[df$pcgroup %in% id.group,]
#     if (sum(as.character(group.df$Add.name) %in% Add.Id)>0) {
#       valid.cmps <- as.character(group.df$Compound[(group.df$Add.name %in% Add.Id)])
#       dx <- group.df[group.df$Compound %in% valid.cmps,]
#       peak.ids <- unique(group.df$Peak.Id[!(group.df$Peak.Id %in% dx$Peak.Id)])
#       return(rbind(dx, group.df[group.df$Peak.Id %in% peak.ids,]))
#     } else {
#       return(group.df)
#     }
#
#   },.parallel = TRUE)
#   return(MH.df)
# }




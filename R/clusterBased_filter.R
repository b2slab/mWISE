#' @name filtering-funs
#' @aliases clusterBased.filter
#' @title Functions to apply cluster-based filtering
#' @description
#' Function \code{clusterBased.filter} selects only those 
#' proposed compounds 
#' for which a quasimolecular adduct or fragment has been also 
#' proposed in another peak of the same cluster.
#' @param df
#' Columns may contain: "Compound", "Add.Id", "Isotope"
#' "Compound" for the proposed candidates.
#' "Add.Id" for the adduct or fragment proposed.
#' "Isotope" to identify the proposed isotopologues.
#' @param Freq
#' Minimum observed frequency to consider an adduct or a fragment 
#' to apply the filter (Def: 0.5).
#' @param Add.Id
#' It indicates the adduct(s) or fragment(s) that are 
#' required to exist.
#' If NULL, those adducts with an observed frequency equal or 
#' higher than 0.50 will be used.
#' @param polarity
#' Acquisition mode of the study. It can be "positive" 
#' or "negative".
#' @param Info.Add
#' Data frame with adducts and in source fragments information. 
#' If NULL, the default mWISE table will be loaded. 
#' The columns should be:
#' \itemize{
#' \item{name} with the name of the adduct or fragment
#' \item{nmol}{ with the number of molecules (i.e., 2M+H: nmol=2 )}
#' \item{charge}{ with the charge of the adduct or 
#' fragment (i.e., M+3H: charge=3)}
#' \item{massdiff}{ with the mass difference 
#' (i.e., M+H: massdiff=1.007276)}
#' \item{quasi}{ with a 1 if the adduct should be considered 
#' as quasi-molecular and a 0, otherwise (Optional)}
#' \item{polarity}{ with a character vector indicating the polarity 
#' of the adduct or fragment. The options are "positive" 
#' or "negative"}}
#' @return
#' Function \code{clusterBased.filter} returns a data frame of 
#' filtered candidates.
#' @examples 
#' data("sample.keggDB")
#' Cpd.Add <- CpdaddPreparation(KeggDB = sample.keggDB, do.Par = FALSE)
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
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply
#' @importFrom utils data

clusterBased.filter <- function(df, Add.Id = NULL, Freq = 0.50, Info.Add = NULL,
                                polarity, do.Par = TRUE, nClust = 2) {
  if (do.Par){
    doParallel::registerDoParallel(nClust)
  }
  if (is.null(Info.Add)){
    data("Info.Add", envir=environment())
  }
  
  if (is.null(Add.Id)){
    Add.Id.freq <- as.character(
      Info.Add$name[(Info.Add$Freq>=Freq)&(Info.Add$polarity %in% polarity)])
    Add.Id.QM <- as.character(
      Info.Add$name[(Info.Add$quasi==1)&(Info.Add$polarity %in% polarity)])
    Add.Id <- unique(c(Add.Id.freq, Add.Id.QM))
  }
  cat(paste("Applying the cluster-based filter using ", 
            paste(Add.Id, collapse = ", "), "...", sep = ""))
  MH.df <- plyr::ldply(unique(df$pcgroup), function(id.group){
    group.df <- df[df$pcgroup %in% id.group,]
    valid.cmps <- as.character(
      group.df$Compound[(group.df$Add.name %in% Add.Id)])
    return(group.df[group.df$Compound %in% valid.cmps,])
  },.parallel = do.Par)
  cat("DONE!","\n")
  return(MH.df)
}



#' @name diffusion-funs
#' @aliases diffusion.input
#' @title Functions to apply diffusion in graphs
#' @description
#' Function \code{diffusion.input} computes the diffusion input 
#' score to perform diffusion in graphs.
#' It uses the functions \code{diffusion.input.Binary} and 
#' \code{diffusion.input.probability}
#' @param df
#' Data frame containing the potential candidates. It is 
#' recommended to use the
#' data frame resulted from mWISE clustered-based filtering. 
#' Columns may contain "Peak.Id" for a peak identifier and
#' "Compound" for a KEGG ID.
#' @param input.type
#' Diffusion input type per compound.
#' "binary" 1 if the compound is proposed.
#' "probability" computes the probability of existence of each 
#' compound.
#' @param background
#' Vector containing a list of KEGG identifiers which will be 
#' set to 0 in the diffusion process. This will have an 
#' effect in the normalization process performed when using 
#' the z score. If NULL, the background will be set to all 
#' the compounds available in df.
#' @param Unique.Annotation
#' Logical (only available when input type="binary"). 
#' If TRUE, the binary diffusion input is computed
#' by only considering those peaks with a unique 
#' annotation (Def: FALSE).
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{diffusion.input} returns a list containing the 
#' diffusion input data frame and a character
#' specifying the diffusion input type.
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
#' Input.diffusion <- diffusion.input(df = MH.Tab,
#'                                   input.type = "probability",
#'                                   Unique.Annotation = FALSE,
#'                                   do.Par = FALSE)
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr llply

diffusion.input<-function(df, input.type = "probability", background = NULL, 
                          Unique.Annotation = FALSE,
                          do.Par = TRUE, nClust = 2){
  df.Input <- df
  if(do.Par)
    doParallel::registerDoParallel(nClust)
  if (Unique.Annotation){
    Unique.ids <- plyr::llply(unique(df$Peak.Id), function(id) {
      if (sum(df$Peak.Id %in% id)==1) return(id)
    },.parallel = do.Par) #%>% unlist()
    Unique.ids <- unlist(Unique.ids)
    df.Input <- df.Input[df.Input$Peak.Id %in% Unique.ids,]
  }

  if (input.type == "binary"){
    Diffusion.Input <- diffusion.input.Binary(df = df.Input, do.Par = do.Par)
  }
  if (input.type == "probability") {
    if (Unique.Annotation){
      Diffusion.Input <- diffusion.input.probability(df = df.Input, 
                                                     do.Par = do.Par)
    } else {
      Diffusion.Input <- diffusion.input.probability(df = df.Input, 
                                                     do.Par = do.Par)
    }
  }
  Whole.Coverage <- data.frame(Compound = unique(df$Compound))
  colnames(Whole.Coverage) <- "Compound"
  Diffusion.Input <- merge(Diffusion.Input, Whole.Coverage, 
                           by = "Compound", all = TRUE)
  Diffusion.Input$Diffusion.Input[is.na(Diffusion.Input$Diffusion.Input)] <- 0

  if(!is.null(background)){
    background <- data.frame(Compound = unique(background))
    Diffusion.Input<-merge(background, Diffusion.Input, 
                           by = "Compound", all = TRUE)
    Diffusion.Input$Diffusion.Input[is.na(Diffusion.Input$Diffusion.Input)]<-0
  }

  return(list(Diffusion.Input=Diffusion.Input,
              input.type= input.type))
}

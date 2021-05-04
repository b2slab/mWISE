#' @name diffusion-funs
#' @aliases set.diffusion
#' @title Functions to apply diffusion in graphs
#' @description
#' Function \code{set.diffusion} diffuses heat in a specific
#' Kernel given a matrix of compounds and its
#' diffusion input score.
#' @param df
#' Object returned by the \code{diffusion.input} function. It is a
#' list containing the diffusion input matrix (data frame 
#' containing a column named "Compound"
#' with the KEGG identifiers of the potential candidates 
#' and a column named "Diffusion.Input") and
#' a character vector indicating the diffusion input type.
#' @param graph
#' Diffusion graph where nodes correspond to KEGG compounds.
#' If NULL, the diffusion graph indicated in \code{graph.name} 
#' will be loaded.
#' @param graph.name
#' Name of the diffusion graphs available in mWISE. 
#' The options are "fella", "RClass3levels" or 
#' "RClass2levels" (Def: "fella").
#' @param K
#' Regularised Laplacian kernel. If NULL, it will be computed 
#' using the \code{regularisedLaplacianKernel} function
#' from DiffuStats R package.
#' @param scores
#' Method of diffusion. Def: c("raw", "ber_s", "z")
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{set.diffusion} returns a list containing 
#' the diffusion results, the compounds discarded during
#' the diffusion process, the compounds present in the network 
#' and the background used.
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
#' data("sample.graph")
#' gMetab <- igraph::as.undirected(sample.graph)
#' diff.Cpd <- set.diffusion(df = Input.diffusion,
#'                           scores = "z",
#'                           graph = gMetab,
#'                           do.Par = FALSE)
#' @export
#' @importFrom igraph as.undirected V
#' @importFrom diffuStats regularisedLaplacianKernel diffuse
#' @importFrom doParallel registerDoParallel
#' @importFrom stats setNames
#' @importFrom utils data

set.diffusion<-function(df, graph = NULL, graph.name = "fella",
                        K = NULL, scores = c("raw", "ber_s", "z"),
                        do.Par = TRUE, nClust = 2){
  if(do.Par) doParallel::registerDoParallel(nClust)
  cat("Preparing diffusion kernel...")
  if (is.null(graph)){
    if (graph.name %in% "fella"){
      Fella <- NULL; data("Fella", envir=environment())
      graph <- igraph::as.undirected(Fella)
    } else if (graph.name %in% "RClass3levels"){
      graph.Rclass <- NULL; data("graph.Rclass", envir=environment())
      graph <- igraph::as.undirected(graph.Rclass)
    } else {
      graph.Rclass.2levels <- NULL
      data("graph.Rclass.2levels", envir=environment())
      graph <- graph.Rclass.2levels
    }}
  if (is.null(K)){
    K <- diffuStats::regularisedLaplacianKernel(graph)
  }
  cat("DONE!","\n"); list.methods <- scores
  cat("Preparing peaks for diffusion:","\n"); meta.params<-df$input.type
  df <- df$Diffusion.Input; Background<-as.character(unique(df$Compound))
  cat(length(Background),"- Compounds as input || " )
  subs.df <- subset(df, !(df$Compound %in% igraph::V(graph)$name))
  Discarded.Compounds <- as.character(unique(subs.df$Compound))
  cat(length(Discarded.Compounds),"- Compounds aren't in the Kernel","\n")
  df<-subset(df, df$Compound %in% igraph::V(graph)$name); df_input<-df
  cat("Implementing diffusion on:",dim(df_input)[1],"peaks...")
  hits <- setNames(df_input$Diffusion.Input, df_input$Compound)
  scores_diff <- lapply(setNames(list.methods, list.methods), function(method){
      diffuStats::diffuse(K = K, scores = hits, method = method)
  })
  scores_diff <- as.data.frame(scores_diff); cat("DONE!","\n")
  cat("Formatting diffusion output...")
  scores_orig <- setNames(rep(NA_real_,nrow(scores_diff)),rownames(scores_diff))
  scores_orig[names(hits)] <- hits; ans <- cbind(scores_orig, scores_diff)
  colnames(ans)[1] <- "Input.Score"
  if(any(is.na(ans[,1]))){
    ans<-ans[-which(is.na(ans[,1])),]
  }
  Result<-data.frame(Compound = rownames(ans));rownames(Result)<-Result$Compound
  Result<- merge(Result,ans, by =0); Result$Row.names<-NULL; cat("DONE!","\n")
  return(list(Diffusion.Results = Result, 
              Discarded.Compounds = Discarded.Compounds,
              Network.Compounds = igraph::V(graph)$name,
              Background = Background))
}

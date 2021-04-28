#' @name diffusion-funs
#' @aliases set.diffusion
#' @title Functions to apply diffusion in graphs
#' @description
#' Function \code{set.diffusion} diffuses heat in a specific
#' Kernel given a matrix of compounds and its
#' diffusion input score.
#' @param df
#' Object returned by the \code{diffusion.input} function. It is a
#' list containing the diffusion input matrix (data frame containing a column named "Compound"
#' with the KEGG identifiers of the potential candidates and a column named "Diffusion.Input") and
#' a character vector indicating the diffusion input type.
#' @param graph
#' Diffusion graph where nodes correspond to KEGG compounds.
#' If NULL, the diffusion graph indicated in \code{graph.name} will be loaded.
#' @param graph.name
#' Name of the diffusion graphs available in mWISE. The options are "fella", "RClass3levels" or "RClass2levels" (Def: "fella").
#' @param K
#' Regularised Laplacian kernel. If NULL, it will be computed using the \code{regularisedLaplacianKernel} function
#' from DiffuStats R package.
#' @param scores
#' Method of diffusion. Def: c("raw", "ber_s", "z")
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{set.diffusion} returns a list containing the diffusion results, the compounds discarded during
#' the diffusion process, the compounds present in the network and the background used.
#' @export

set.diffusion<-function(df, graph = NULL, graph.name = "fella",
                        K = NULL, scores = c("raw", "ber_s", "z"),
                        do.Par = TRUE, nClust = 2){

  if(do.Par)
    doParallel::registerDoParallel(nClust)

  cat("Preparing diffusion kernel...")
  # Load graph (network)
  if (is.null(graph)){
    #data(graph)
    if (graph.name %in% "fella"){
      graph <- igraph::as.undirected(mWISE::Fella)
    } else if (graph.name %in% "RClass3levels"){
      graph <- igraph::as.undirected(mWISE::graph.Rclass)
    } else {
      graph <- mWISE::graph.Rclass.2levels
    }

  }
  g.metab<-graph

  # Compute graph kernel
  if (is.null(K)){
    K <- diffuStats::regularisedLaplacianKernel(g.metab)
  }
  cat("DONE!","\n")

  # Scores to try
  list.methods <- scores

  cat("Preparing peaks for diffusion:","\n")
  # Filter Compounds that are in the network
  meta.params<-df$input.type
  df <- df$Diffusion.Input
  Background<-as.character(unique(df$Compound))
  cat(length(Background),"- Compounds as input || " )
  # Discarded.Compounds<-df %>%
  #   subset(!Compound %in% V(g.metab)$name) %>%
  #   .$Compound %>% unique() %>% as.character()
  subs.df <- subset(df, !(df$Compound %in% igraph::V(g.metab)$name))
  Discarded.Compounds <- as.character(unique(subs.df$Compound))
  cat(length(Discarded.Compounds),"- Compounds aren't in the Kernel","\n")
  #df<-df %>% subset(Compound %in% V(g.metab)$name)
  df<-subset(df, df$Compound %in% igraph::V(g.metab)$name)
  df_input<-df
  cat("Implementing diffusion on:",dim(df_input)[1],"peaks...")
  # if (any(meta.params == c("binary"))){
  #   cat("\n","WARNING! ",meta.params,"input is not suited for z, ber_p and mc diffusion scores..." ,"\n")
  #   df_input$Diffusion.Input <-df_input$Diffusion.Input +
  #     rnorm(length(df_input$Diffusion.Input), mean =0, sd = 1e-6) %>% abs()
  #   cat("Gaussian noise added to input. Output may be missrepresented for: z, ber_p and mc scores" ,"\n")
  # }

  # Prepare diffusion - the input must be a named vector
  # (assuming the names are in the network)
  hits <- setNames(df_input$Diffusion.Input, df_input$Compound)
  # Diffusion process (one for each desired score)
  scores_diff <- sapply(
    setNames(list.methods, list.methods),
    function(method) {
      diffuStats::diffuse(K = K, scores = hits, method = method)
    }
  )
  cat("DONE!","\n")

  cat("Formatting diffusion output...")
  # Orininal scores (which are the number of hits)
  scores_orig <- setNames(rep(NA_real_, nrow(scores_diff)), rownames(scores_diff))
  scores_orig[names(hits)] <- hits

  # Place everything in a matrix to export
  ans <- cbind(scores_orig, scores_diff)
  colnames(ans)[1] <- "Input.Score"
  if(any(is.na(ans[,1]))){
    ans<-ans[-which(is.na(ans[,1])),]
  }
  Result<-data.frame(Compound = rownames(ans))
  rownames(Result)<-Result$Compound
  Result<- merge(Result,ans, by =0); Result$Row.names<-NULL
  cat("DONE!","\n")

  return(list(Diffusion.Results = Result,
              Discarded.Compounds = Discarded.Compounds,
              Network.Compounds =igraph::V(g.metab)$name,
              Background = Background))
}

#' @name performanceEvaluation
#' @title Function to evaluate mWISE performance
#' @description
#' Function \code{performanceEvaluation} computes performance 
#' metrics when a reference table is available.
#' @param Annotated.dataset
#' Object returned by \code{mWISE.annotation} function.
#' @param df.Ref
#' Reference data frame. It should contain a column named "Peak.Id"
#' with a peak identifier that corresponds to the mWISE Peak 
#' identifier, and a column named "Kegg"
#' with the reference KEGG identifiers.
#' @param top.cmps
#' Number of top compounds to be considered in performance 
#' evaluation.
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE.
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{performanceEvaluation} returns a data frame 
#' with the performance of each mWISE stage.
#' @examples 
#' data("sample.keggDB")
#' Cpd.Add <- CpdaddPreparation(KeggDB = sample.keggDB, do.Par = FALSE)
#' data(sample.dataset)
#' Peak.List <- sample.dataset$Negative$Input
#' Intensity.idx <- seq(27,38)
#' data("sample.graph")
#' gMetab <- igraph::as.undirected(sample.graph)
#' Annotated.List <- mWISE.annotation(Peak.List = Peak.List,
#'                                    polarity = "negative",
#'                                    diffusion.input.type = "binary",
#'                                    score = "raw",
#'                                    Cpd.Add = Cpd.Add,
#'                                    graph = gMetab,
#'                                    Unique.Annotation = TRUE,
#'                                    Intensity.idx = Intensity.idx,
#'                                    do.Par = FALSE)
#'  df.Ref <- sample.dataset$Negative$Output
#'  performanceEvaluation(Annotated.dataset = Annotated.List, 
#'                        df.Ref = df.Ref, top.cmps = 3)
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply

performanceEvaluation <- function(Annotated.dataset, df.Ref, top.cmps, 
                                  do.Par = TRUE, nClust){
  if (do.Par) doParallel::registerDoParallel(nClust)
  Peak.Cpd <- Annotated.dataset$Annotated.Tab
  merge.df <- merge(df.Ref[,c("Peak.Id", "Kegg")],
                    Peak.Cpd, by = "Peak.Id")
  stages <- c("Matching", "Filtering", "Diffusion")
  Annotated.dataset$Ranked.Tab <- Annotated.dataset$Ranked.Tab[
    Annotated.dataset$Ranked.Tab$Ranking<=top.cmps,]
  metrics.stages <- plyr::ldply(seq_len(3), function(j){
    i <- c(1,3,5)[j]; df <- Annotated.dataset[[i]]
    not.Annotated <- sum(!(unique(df.Ref$Peak.Id) %in% unique(df$Peak.Id)))
    Metrics <- plyr::ddply(merge.df, "Peak.Id", function(dy){
      dx <- df[df$Peak.Id %in% unique(dy$Peak.Id),]
      Ref.Cmp <- unique(as.character(dy$Kegg))
      if (nrow(dx)>0){
        TP <- sum(dx$Compound %in% Ref.Cmp); TP <- ifelse(TP>0, 1, 0)
        FP <- sum(!(dx$Compound %in% Ref.Cmp))
        FN <- sum(!(Ref.Cmp %in% dx$Compound))
        TN <- sum(
          (!(dy$Compound %in% dx$Compound)) & (!(dy$Compound %in% Ref.Cmp)))
        data.frame(TP = TP, FP = FP, TN = TN, FN = FN, 
                   Background = sum(TP, FP, TN, FN))
      } else {
        TP <- 0; FP <- 0; FN <- 1
        Background <- nrow(dy) + ifelse(sum(dy$Compound %in% Ref.Cmp)==0,1,0)
        TN <- Background-FN
        data.frame(TP = TP, FP = FP, TN = TN, FN = FN, 
                   Background = sum(TP, FP, TN, FN))
      }
    },.parallel = do.Par)
    Total.Metrics <- as.data.frame(
      t(colSums(Metrics[,c("TP", "FP", "TN", "FN")])))
    Total.Metrics$NotAnnotated <- not.Annotated
    Total.Metrics$Ref.Numbers <- length(unique(df.Ref$Peak.Id))
    Total.Metrics$Sensitivity <-
      round(Total.Metrics$TP/(Total.Metrics$TP+Total.Metrics$FN), 2)
    Total.Metrics$Specificity <-
      round(Total.Metrics$TN/(Total.Metrics$TN+Total.Metrics$FP), 2)
    Total.Metrics$Precision <-
      round(Total.Metrics$TP/(Total.Metrics$TP+Total.Metrics$FP), 2)
    Total.Metrics$Accuracy <-
      round((Total.Metrics$TP+Total.Metrics$TN)/
      (Total.Metrics$TP+Total.Metrics$TN+Total.Metrics$FP+Total.Metrics$FN), 2)
    Total.Metrics$F1 <- round(Total.Metrics$TP/(Total.Metrics$TP+(0.5*(
        Total.Metrics$FP+Total.Metrics$FN))),2)
    Total.Metrics <- data.frame(Stage = stages[j], Total.Metrics)
    return(Total.Metrics)
  }); return(metrics.stages)
}

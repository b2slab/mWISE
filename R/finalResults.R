#' @name diffusion-funs
#' @aliases finalResults
#' @aliases modifiedTabs
#' @title Functions to apply diffusion in graphs
#' @description
#' Function \code{finalResults} prepares the final table ranked by the diffusion scores computed.
#' Function \code{modifiedTabs} prepares the tables that result from the matching stage and the
#' filtering stage to evaluate their performance.
#' @example
#' Ranked.Tab <- finalResults(Diff.Tab = Diff.Tab, score = "z",
#'                            do.Par = TRUE, nClust = detectCores()-1)
#' @param Diff.Tab
#' Data frame that results from the diffusion step.
#' @param score
#' Method of diffusion. Possible values are: "raw", "ber_s" and "z"
#' @param df
#' Data frame containing the results from the matching or the filtering stage.
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE.
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{finalResults} returns a table containing the final potential candidates ranked by diffusion score.
#' Function \code{modifiedTabs} returns a table containing the potential candidates in a specific mWISE stage without repeated peaks.
#' @export

finalResults <- function(Diff.Tab, score, do.Par = TRUE, nClust){
  doParallel::registerDoParallel(nClust)
  cat("Preparing the diffusion-based ranked table...")
  Diff.Tab <- plyr::ldply(unique(Diff.Tab$Peak.Id), function(id){
    dx <- Diff.Tab[Diff.Tab$Peak.Id %in% id,]
    if (sum(duplicated(dx$Compound))>0){
      cmps.rep <- as.character(dx$Compound[duplicated(dx$Compound)])
      dn <- dx[!(dx$Compound %in% cmps.rep),]
      dr <- plyr::ldply(cmps.rep, function(c){
        dy <- dx[dx$Compound %in% c,]
        Add <- paste(dy$Add.name, collapse = " | ")
        dy <- dy[1,]
        dy$Add.name <- Add
        return(dy)
      })
      dx <- rbind(dn,dr)
      return(dx)
    } else {
      return(dx)
    }
  },.parallel = do.Par)

  Ranked.Tab <- plyr::ldply(unique(Diff.Tab$Peak.Id), function(id){
    dx <- Diff.Tab[Diff.Tab$Peak.Id %in% id,]
    if (nrow(dx)>1){
      if (sum(is.na(dx[,score])) < nrow(dx)){
        dz <- dx[!is.na(dx[,score]),]
        dz$Ranking <- rank(-dz[,score])
        return(dz)
      } else {
        dx$Ranking <- nrow(dx)
        return(dx)
      }
    } else {
      dx$Ranking <- 1
      return(dx)
    }
  },.parallel = do.Par)
  cat("DONE!","\n")
  return(Ranked.Tab)
}

modifiedTabs <- function(df, do.Par = TRUE, nClust){
  df <- plyr::ldply(unique(df$Peak.Id), function(id){
    dx <- df[df$Peak.Id %in% id,]
    if (sum(duplicated(dx$Compound))>0){
      cmps.rep <- as.character(dx$Compound[duplicated(dx$Compound)])
      dn <- dx[!(dx$Compound %in% cmps.rep),]
      dr <- plyr::ldply(cmps.rep, function(c){
        dy <- dx[dx$Compound %in% c,]
        Add <- paste(dy$Add.name, collapse = " | ")
        dy <- dy[1,]
        dy$Add.name <- Add
        return(dy)
      })
      dx <- rbind(dn,dr)
      return(dx)
    } else {
      return(dx)
    }
  },.parallel = TRUE)
}












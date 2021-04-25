evaluation.performance <- function(Original.Matching, df, df.Ref,
                                   do.Par){
  df.Ref.peaks <- df.Ref$Peak.Id %>% unique() %>% length()
  Total.df <- merge(df.Ref[,c("Kegg", "Peak.Id")],
                    df, by = c("Peak.Id"))
  Original.df <- merge(df.Ref[,c("Kegg", "Peak.Id")],
                       Original.Matching, by = "Peak.Id")

  not.Annotated <- sum(!(unique(df.Ref$Peak.Id) %in% unique(Total.df$Peak.Id)))
  Metrics <- ddply(Total.df, "Peak.Id", function(dx){
    dy <- Original.df[Original.df$Peak.Id %in% unique(dx$Peak.Id),]
    Ref.Cmp <- unique(as.character(dx$Kegg))
    TP <- (dx$Compound %in% Ref.Cmp) %>% sum()
    FP <- (!(dx$Compound %in% Ref.Cmp)) %>% sum()
    FN <- (!(Ref.Cmp %in% dx$Compound)) %>% sum()
    TN <- ((!(dy$Compound %in% dx$Compound)) & (!(dy$Compound %in% Ref.Cmp))) %>% sum()
    data.frame(TP = TP, FP = FP, TN = TN, FN = FN, Background = nrow(dy))
  },.parallel = do.Par)
  Total.Metrics <- colSums(Metrics[,c("TP", "FP", "TN", "FN")]) %>% t() %>% as.data.frame()
  Total.Metrics$NotAnnotated <- not.Annotated
  Total.Metrics$Ref.Peaks <- df.Ref.peaks
  Total.Metrics$Sensitivity <- Total.Metrics$TP/(Total.Metrics$TP+Total.Metrics$FN+Total.Metrics$NotAnnotated)
  Total.Metrics$Specificity <- Total.Metrics$TN/(Total.Metrics$TN+Total.Metrics$FP)
  Total.Metrics$Precision <- Total.Metrics$TP/(Total.Metrics$TP+Total.Metrics$FP)
  Total.Metrics$Accuracy <- (Total.Metrics$TP+Total.Metrics$TN)/(Total.Metrics$TP+Total.Metrics$TN+Total.Metrics$FP+Total.Metrics$FN)

  return(Total.Metrics)
}

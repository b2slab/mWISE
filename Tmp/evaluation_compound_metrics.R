evaluation.compound.metrics <- function(Original.Matching, Peak.Cpd, df.Ref) {
  Ref.Cmps <- df.Ref$Kegg %>% unique() %>% as.character()
  Total.Refs <- length(Ref.Cmps)
  TP <- sum(as.character(unique(Peak.Cpd$Compound)) %in% Ref.Cmps)
  FP <- sum(!(as.character(unique(Peak.Cpd$Compound)) %in% Ref.Cmps))
  FN <- sum(!(Ref.Cmps %in% as.character(unique(Peak.Cpd$Compound))))
  TN <- sum((!(as.character(unique(Original.Matching$Compound)) %in% Ref.Cmps)) &
    (!as.character(unique(Original.Matching$Compound)) %in% as.character(unique(Peak.Cpd$Compound))))

  Total.Metrics <- data.frame(TP = TP, FP = FP, TN = TN, FN = FN)
  Total.Metrics$Ref.Compounds <- Total.Refs
  Total.Metrics$Sensitivity <- Total.Metrics$TP/(Total.Metrics$TP+Total.Metrics$FN)
  Total.Metrics$Specificity <- Total.Metrics$TN/(Total.Metrics$TN+Total.Metrics$FP)
  Total.Metrics$Precision <- Total.Metrics$TP/(Total.Metrics$TP+Total.Metrics$FP)
  Total.Metrics$Accuracy <- (Total.Metrics$TP+Total.Metrics$TN)/(Total.Metrics$TP+Total.Metrics$TN+Total.Metrics$FP+Total.Metrics$FN)
  return(Total.Metrics)
}

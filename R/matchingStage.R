#' @name matching-funs
#' @aliases matchingStage
#' @title Functions to match a peak-intensity table to KEGG database
#' @description
#' Function \code{matchingStage} performs a fast matching of the peak-intensity table to KEGG database
#' using adducts and fragments knowledge.
#' @examples
#' matchingStage <- function(Peak.List, polarity = "positive", 
#' do.Par = TRUE, nClust = 2)
#' @param Peak.List
#' Data frame containing the LC-MS features.
#' Columns should contain:
#' \itemize{
#' \item{Peak.Id}{ for a peak identifier}
#' \item{mz}{ for a mass-to-charge ratio value}
#' \item{rt}{ for the retention time}
#' \item{Intensities for each sample}
#' }
#' @param force.mass.range
#' Logical. If TRUE, the m/z range is computed using the user-defined ppm tolerance or mz.range (Def: TRUE).
#' @param mass.range.type
#' A character indicating if the m/z range computation is performed using tolerance in ppm ("ppm.mode").
#' or m/z uncertainty ("mz.range") (Def: "ppm.mode").
#' @param mz.range
#' If mass.range.type = "mz.range", it indicates the m/z range uncertainty (i.e. mz.range = 0.001).
#' @param ppm
#' If mass.range.type = "ppm.mode", it indicates the tolerance in ppm (i.e. ppm = 10).
#' @param polarity
#' Acquisition mode of the study. It can be "positive" or "negative".
#' @param Add.List
#' List of adducts or fragments to consider.
#' @param Cpd.Add
#' Compound to adduct matrix. It can be built using \code{CpdaddPreparation} function.
#' If NULL, the default information will be used.
#' @param do.Par
#' Logical. If parallel computing is required (Def: TRUE).
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{matchingStage} returns a list containing the input peak list and the table of annotated peaks,
#' containing all the potential KEGG candidates for each LC-MS feature.
#' @export

matchingStage <- function(Peak.List, force.mass.range = TRUE, mass.range.type = "ppm.mode",
                          mz.range = NULL, ppm = 10, polarity = "positive",
                          Add.List = NULL, Cpd.Add,
                          do.Par = TRUE, nClust){

  # Formatting Peak.List
  cat("Matching the input peak list to KEGG database...")
  Peak.List <- mass.range.format(Peak.List = Peak.List, force.mass.range = force.mass.range,
                                 mass.range.type = mass.range.type, mz.range = mz.range, ppm = ppm)

  # if (is.null(Cpd.Add) & is.null(Add.List)) {
  #   Cpd.Add <- mWISE::Cpd.Add
  # } else if (is.null(Cpd.Add) & (!is.null(Add.List))) {
  #   Cpd.Add <- mWISE::Cpd.Add
  #   Cpd.Add <- Cpd.Add[as.character(Cpd.Add$Add.name) %in% Add.List,]
  # } else{
  #   Cpd.Add <- Cpd.Add
  # }
  if (!is.null(Add.List)){
    Cpd.Add <- Cpd.Add[as.character(Cpd.Add$Add.name) %in% Add.List,]
  }
  
  Cpd.Add <- Cpd.Add[Cpd.Add$Polarity %in% polarity,]

  if (max(Peak.List$mzmax)>2000){
    # m/z<2000
    Peak.List1 <- Peak.List[Peak.List$mzmax<=2000,]
    Cpd.Add1 <- subset(Cpd.Add, Cpd.Add$mz.Add<=2000)
    gr.PeakList1 <- GenomicRanges::GRanges(seqnames = 'dummy', strand = '*',
                            ranges = IRanges::IRanges(start = round(Peak.List1$mzmin*10^6, 0),
                                             end = round(Peak.List1$mzmax*10^6, 0)))
    gr.CpdAdd1 <- GenomicRanges::GRanges(seqnames = 'dummy', strand = "*",
                          ranges = IRanges::IRanges(start = round(Cpd.Add1$mz.Add*10^6, 0), width = 1))
    overlaps1<-data.frame(GenomicRanges::findOverlaps(query = gr.PeakList1, subject = gr.CpdAdd1))
    Peak.Cpd1 <- cbind(Peak.List1[overlaps1$queryHits,c("Peak.Id", "mz", "rt")],
                      Cpd.Add1[overlaps1$subjectHits, c("Compound", "exact_mass", "Add.name")])
    # m/z>2000
    Peak.List2 <- Peak.List[Peak.List$mz>2000,]
    Cpd.Add2 <- subset(Cpd.Add, Cpd.Add$mz.Add>2000)
    gr.PeakList2 <- GenomicRanges::GRanges(seqnames = 'dummy', strand = '*',
                            ranges = IRanges::IRanges(start = round(Peak.List2$mzmin*10^4, 0),
                                             end = round(Peak.List2$mzmax*10^4, 0)))
    gr.CpdAdd2 <- GenomicRanges::GRanges(seqnames = 'dummy', strand = "*",
                          ranges = IRanges::IRanges(start = round(Cpd.Add2$mz.Add*10^4, 0), width = 1))
    overlaps2<-data.frame(GenomicRanges::findOverlaps(query = gr.PeakList2, subject = gr.CpdAdd2))

    Peak.Cpd2 <- cbind(Peak.List2[overlaps2$queryHits,c("Peak.Id", "mz", "rt")],
                      Cpd.Add2[overlaps2$subjectHits, c("Compound", "exact_mass", "Add.name")])
    Peak.Cpd <- rbind(Peak.Cpd1, Peak.Cpd2)

  } else {
    # Formatting Peak.List
    Cpd.Add <- subset(Cpd.Add, Cpd.Add$mz.Add<=2000)
    gr.PeakList <- GenomicRanges::GRanges(seqnames = 'dummy', strand = '*',
                           ranges = IRanges::IRanges(start = round(Peak.List$mzmin*10^6, 0),
                                            end = round(Peak.List$mzmax*10^6, 0)))
    gr.CpdAdd <- GenomicRanges::GRanges(seqnames = 'dummy', strand = "*",
                         ranges = IRanges::IRanges(start = round(Cpd.Add$mz.Add*10^6, 0), width = 1))
    overlaps<-data.frame(GenomicRanges::findOverlaps(query = gr.PeakList, subject = gr.CpdAdd))

    Peak.Cpd <- cbind(Peak.List[overlaps$queryHits,c("Peak.Id", "mz", "rt")],
                      Cpd.Add[overlaps$subjectHits, c("Compound", "exact_mass", "Add.name")])
  }
  cat("DONE!","\n")

  return(list(Peak.List = Peak.List,
              Peak.Cpd = Peak.Cpd))
}

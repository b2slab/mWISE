#' @name mWISE.annotation
#' @title Function to perform the complete annotation pipeline of mWISE
#' @description
#' Wrapper function that performs the complete workflow of mWISE to annotate a peak-intensity matrix.
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
#' @param Intensity.idx
#' Numeric vector indicating the column index for the intensities
#' @param Rt.05
#' Retention time value to get a similarity of 0.5.
#' @param Freq
#' Minimum observed frequency to consider an adduct or a fragment to apply the cluster-based filter.
#' @param Add.Id
#' It indicates the adduct(s) or fragment(s) that are considered to apply the cluster-based filter.
#' If NULL, those adducts with an observed frequency equal or higher than 0.10 will be used.
#' @param diffusion.input.type
#' Diffusion input type per compound.
#' "binary" 1 if the compound is proposed.
#' "probability" computes the probability of existence of each compound.
#' @param background
#' Vector containing a list of KEGG identifiers which will be set to 0 in the diffusion process. This will have an
#' effect in the normalization process performed when using the z score.
#' If NULL, the background will be set to all the compounds available in df.
#' @param Unique.Annotation
#' Logical (only available when input type="binary"). If TRUE, the binary diffusion input is computed
#' by only considering those peaks with a unique annotation (Def: FALSE).
#' @param graph
#' Diffusion graph where nodes correspond to KEGG compounds.
#' If NULL, the diffusion graph indicated in \code{graph.name} will be loaded.
#' @param graph.name
#' Name of the diffusion graphs available in mWISE. The options are "fella", "RClass3levels" or "RClass2levels" (Def: "fella").
#' @param K
#' Regularised Laplacian kernel. If NULL, it will be computed using the \code{regularisedLaplacianKernel} function
#' from DiffuStats R package.
#' @param score
#' Method of diffusion. Def: c("raw", "ber_s", "z")
#' @param do.Par
#' TRUE if parallel computing is required. Def: TRUE
#' @param nClust
#' Number of clusters that may be used. Def: Number of clusters - 1.
#' @return
#' Function \code{matchingStage} returns a list containing the input peak list and the table of annotated peaks,
#' containing all the potential KEGG candidates for each LC-MS feature.
#' @export

mWISE.annotation <- function(Peak.List, force.mass.range = TRUE, mass.range.type = "ppm.mode",
                             mz.range = NULL, ppm = 10, polarity = "positive", Add.List = NULL, Cpd.Add,
                             Intensity.idx, Rt.05 = 5, Freq = 0.50, Add.Id = NULL,
                             background = NULL,  diffusion.input.type = "probability", Unique.Annotation = FALSE,
                             graph = NULL, K = NULL, score = "z", graph.name = "fella",
                             nClust = 2, do.Par = TRUE) {

  if (do.Par){
    doParallel::registerDoParallel(nClust)
  }
  
  Annotated.List <- matchingStage(Peak.List = Peak.List, force.mass.range = force.mass.range,
                                  mass.range.type = mass.range.type, mz.range = mz.range,
                                  ppm = ppm, polarity = polarity, Add.List = Add.List,
                                  Cpd.Add = Cpd.Add, do.Par = do.Par, nClust = nClust)

  Annotated.Tab <- Annotated.List$Peak.Cpd

  clustered <- featuresClustering(Peak.List = Peak.List,
                                  Intensity.idx = Intensity.idx, Rt.05 = Rt.05,
                                  do.Par = do.Par, nClust = nClust)
  Annotated.Tab <- merge(Annotated.Tab,
                         clustered$Peak.List[,c("Peak.Id", "pcgroup")],
                         by = "Peak.Id")
  # M+H filter
  # if (is.null(Info.Add)){
  #   Info.Add <- mWISE::Info.Add
  # }

  #QM.Adducts <- as.character(Info.Add$name[(Info.Add$quasi==1) & (Info.Add$polarity %in% polarity)])
  MH.Tab <- clusterBased.filter(df = Annotated.Tab, Add.Id = Add.Id, Freq = Freq, polarity = polarity)

  Annotated.Tab <- modifiedTabs(df = Annotated.Tab, do.Par = do.Par,
                                nClust = nClust)
  MH.Tab <- modifiedTabs(df = MH.Tab, do.Par = do.Par, nClust = nClust)

  # Diffusion Input and set diffusion
  Input.diffusion <- diffusion.input(df = MH.Tab, background = background,
                                     input.type = diffusion.input.type,
                                     Unique.Annotation = Unique.Annotation,
                                     do.Par = do.Par, nClust = nClust)


  diff.Cpd <- set.diffusion(df = Input.diffusion,
                            scores = score, graph = graph, K = K,
                            graph.name = graph.name,
                            do.Par = do.Par, nClust = nClust)
  Diffusion.Results <- diff.Cpd$Diffusion.Results
  MH.Tab <- recoveringPeaks(Annotated.Tab = Annotated.Tab, MH.Tab = MH.Tab)
  Diff.Tab <- merge(x = MH.Tab, y = Diffusion.Results,
                    by = "Compound", all.x = TRUE)

  Ranked.Tab <- finalResults(Diff.Tab = Diff.Tab, score = score)

  return(list(Annotated.Tab = Annotated.Tab,
              Clustered.Tab = clustered,
              MH.Tab = MH.Tab,
              Diff.Tab = Diff.Tab,
              Ranked.Tab = Ranked.Tab))

}

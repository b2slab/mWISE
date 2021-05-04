#' @name matching-funs
#' @aliases mass.range.format
#' @title Functions to match a peak-intensity table to 
#' KEGG database
#' @description
#' Function \code{mass.range.format} builds mz ranges 
#' for a peak-intensity table using \code{mz.ppm.converter} 
#' if required.
#' @param Peak.List
#' Columns may be: "peak.id","mz","rt","mzmin","mzmax"
#' "peak.id" for a peak identifier
#' "mz" for detected mass.
#' "mzmin" for minimumum of the range of mass detected.
#' "mzmax" for maximum of the range of mass detected.
#' "rt" for the eluting time.
#' @param force.mass.range
#' Logical. If TRUE the mass range is computed 
#' using the user-defined ppm tolerance or 
#' mz.range (Def: TRUE).
#' @param mass.range.type
#' A character indicating if the m/z range 
#' computation is performed
#' using tolerance in ppm ("ppm.mode") or 
#' m/z uncertainty ("mz.range").
#' Def: "ppm.mode".
#' @param mz.range
#' If mass.range.type = "mz.range", it indicates 
#' the m/z range uncertainty,
#' (i.e. mz.range = 0.001).
#' @param ppm
#' If mass.range.type = "ppm.mode", it 
#' indicated the tolerance 
#' in ppm, i.e. ppm = 10.
#' @return
#' Function \code{mass.range.format} returns 
#' the input peak list 
#' with two additional columns containing the 
#' mzmin and mzmax values.

mass.range.format <- function(Peak.List,
                              force.mass.range = TRUE,
                              mass.range.type = "ppm.mode",
                              mz.range, ppm = 10){
  # Formating Peak List.
  if (!(any(colnames(Peak.List) =="mzmin") & 
        any(colnames(Peak.List) =="mzmax")) | force.mass.range){
    Peak.List<- data.frame(Peak.Id = Peak.List$Peak.Id,
                           mz= Peak.List$mz,
                           rt =Peak.List$rt)

    if (mass.range.type == "mz.range"){
      Peak.List$mzmin<-Peak.List$mz-mz.range
      Peak.List$mzmax<-Peak.List$mz+mz.range
    }

    if (mass.range.type == "ppm.mode"){
      mzrange <- mz.ppm.converter(tol=ppm,  mz=Peak.List$mz, 
                                  conversion=c("ppm.to.mz"))
      Peak.List$mzmin<-Peak.List$mz-mzrange
      Peak.List$mzmax<-Peak.List$mz+mzrange
    }

  } else {
    Peak.List<- data.frame(Peak.Id = Peak.List$Peak.Id,
                           mz= Peak.List$mz,
                           mzmin = Peak.List$mzmin,
                           mzmax = Peak.List$mzmax,
                           rt = Peak.List$rt)
    mz.range <- "experimental"
  }
  return(Peak.List)
}

#' Extract HDM values for relevant cells and resolution for use in RA calculations
#' @param myHDMSpp_NO vector of TAXON_IDs for species to be included in output
#' @param myCropRasters list of rasters and indices and cell values created by function cropNAborder()
#' @param RasterRes numeric raster resolution of the analysis in metres (225 or 75
#'   usually set in settings file or shiny app)
#' @export
makeHDMVals <- function(myHDMSpp_NO = HDMSpp_NO,
                        myCropRasters = cropRasters,
                        RasterRes = myFHAnalysis$RasterRes){
  load(paste0("./HDMS/HDMVals", RasterRes, ".rdata"))
  myHDMVals <- HDMVals[myCropRasters$IDX, as.character(myHDMSpp_NO)]
  return(myHDMVals)
}

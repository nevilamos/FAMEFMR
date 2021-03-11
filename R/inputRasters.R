#' Set correct input general rasters
#'
#' @param RasterRes numeric raster resolution of the analysis in metres (
#'   usually set in settings file or shiny app)
#'
#' @return list of input raster names correct for RasterRes or error if
#'   RasterRes is not 75 or 225
#' @export
inputRasters <- function(RasterRes){
  #General Input Rasters change name depending on RasterRes
  if(RasterRes%in%c(75,225)){
    if (RasterRes == 225){
      y <- list(REGION.tif = "LF_REGION_225.tif",
                EFG.tif = "EFG_NUM_225.tif",
                PLM_GEN.tif = "PLM_GEN_225.tif",
                IDX.tif = "IndexVals225.tif",
                FIREFMZ.tif = "FIRE_FMZ_225.tif",
                DELWP.tif = "DELWP_REGION_225.tif"
      )

    }
    else{
      y <- list(REGION.tif = "LF_REGION_75.tif",
                EFG.tif = "EFG_NUM_75.tif",
                PLM_GEN.tif = "PLM_GEN_75.tif",
                IDX.tif = "IndexVals75.tif",
                FIREFMZ.tif = "FIRE_FMZ_75.tif",
                DELWP.tif = "DELWP_REGION_75.tif"
      )
    }
  }
  else {
    stop("75 and 225 are the only permitted RasterRes values for analysis")
  }
  return(y)
}

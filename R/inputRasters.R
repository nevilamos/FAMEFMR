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
      y <- list(REGION = "LF_REGION_225.tif",
                EFG = "EFG_NUM_225.tif",
                PLM = "PLM_GEN_225.tif",
                IDX = "IndexVals225.tif",
                FIREFMZ = "FIRE_FMZ_225.tif",
                DELWP = "DELWP_REGION_225.tif"
      )

    }
    else{
      y <- list(REGION = "LF_REGION_75.tif",
                EFG = "EFG_NUM_75.tif",
                PLM = "PLM_GEN_75.tif",
                IDX = "IndexVals75.tif",
                FIREFMZ = "FIRE_FMZ_75.tif",
                DELWP = "DELWP_REGION_75.tif"
      )
    }
  }
  else {
    stop("75 and 225 are the only permitted RasterRes values for analysis")
  }
  return(y)
}

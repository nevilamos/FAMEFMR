
#' Calculate all combinations of input raster values
#'
#' @param myFHAnalysis list containing all the fire history spatial attributes created by function fhProcess()
#' @param myCropRasters list of rasters and indices and cell values created by function cropToOutput()
#' @param myRasterRes Resolution of Rasters to be used in the analysis set in settings file in R-script version needs to be set in shiny version
#' @param puPoly either NULL if no planning unit polygons are to be included in the combinations or path to planning unit polygon shapefile covering the entire extent in VicGrid94 projection
#'
#' @return list of:
#' \itemize{
#' \item data.table giving all combinations of cell values from
#'   the input rasters for the FAME analysis
#' \item integer index mapping unique combinations (above) to raster cells
#' }
#' @importFrom  data.table as.data.table setkey setDT
#' @export
calcU_All_Combs<-function (myFHAnalysis = FHAnalysis,
                           myCropRasters = cropRasters,
                           myRasterRes = RasterRes,
                           puPoly = NULL
                           )
						   {

  # get input data
  TimeRange <- as.integer(myFHAnalysis$TimeSpan)
  TimeNames <- as.character(myFHAnalysis$TimeSpan)
  LTR <- length(TimeRange)
  r <- eval(myCropRasters$rasterDef)
  FH_ID <- myCropRasters$FH_ID#terra::values(r)
  gc()

  # get combination input data
  PLM <- myCropRasters$PLM
  EFG <- myCropRasters$EFG
  FIRE_REG <- as.integer(myCropRasters$RGN)
  FIREFMZ <- as.integer(myCropRasters$FIREFMZ)
  DELWP <- myCropRasters$DELWP

  if(is.null(puPoly) == FALSE){
    PU<-terra::values(terra::rasterize(x = terra::vect(puPoly),y = r,field = "PU"))
    # combine into data table
    AllCombs<-data.table::as.data.table(cbind(FH_ID,EFG,FIRE_REG,FIREFMZ,PLM,DELWP,PU))
    rm(PU)
  }else{
      # combine into data table
  AllCombs<-data.table::as.data.table(cbind(FH_ID,EFG,FIRE_REG,FIREFMZ,PLM,DELWP))
    }
  rm(PLM,EFG,FIRE_REG,FIREFMZ,FH_ID,DELWP)
  # clean memory
  gc()

  # calculate the unique combinations of fire history (FH_ID)
  # and EFG and various other admin datasets
  # This then allows the reduction of the subsequent calculation matrices
  # to the number of unique cases of EFG and fire history rather than the number of raster cells
  U_AllCombs <- mgcv::uniquecombs(AllCombs)
  #index of the unique combs for linking back to raster
  Index_AllCombs <- attributes(U_AllCombs)$index
  # remove variables no longer needed
  rm(AllCombs)
  # clean memory
  gc()


  # get the number of pixels in each unique combination
  # this can be used at the end of the process to calculate area cases
  nPixel <- as.vector(table(Index_AllCombs))
  #Add index column a and count of pixels to unique combs matrix
  Index <- 1:nrow(U_AllCombs)
  # create dataframe of all combinations and the pixel count
  U_AllCombs <- cbind(Index, nPixel, U_AllCombs)


  # using data.table formatted left join
  # (the "left" table is the one in brackets on the right)
  data.table::setDT(TFI_LUT)
  data.table::setkey(TFI_LUT, "EFG")
  U_AllCombs <- data.table::as.data.table(U_AllCombs)
  data.table::setkey(U_AllCombs, "EFG")
  U_AllCombs_TFI <- TFI_LUT[U_AllCombs]
  #have to reset the index to Index to return to the original order which is needed for cbinds below.
  data.table::setkey(U_AllCombs_TFI, "Index")
  U_AllCombs_TFI <- Join_Names(U_AllCombs_TFI)
  U_AllCombs_TFI$PLM[U_AllCombs_TFI$PLM == 1] <- "Public Land"
  U_AllCombs_TFI$Hectares <- U_AllCombs_TFI$nPixel * cellsToHectares(RasterMetres = myRasterRes)
  # return function data
  return(list("U_AllCombs_TFI" = U_AllCombs_TFI,"Index_AllCombs" = Index_AllCombs))
}

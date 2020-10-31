
#' Calculate all combinations of input raster values
#'
#' @param FHAnalysis list containing all the fire history spatial attributes created by function fhProcess()
#' @param cropRasters list of rasters and indices and cell values created by function cropNAborder()
#'
#' @return list of:
#'\
#' @export
calcU_All_Combs<- function(FHAnalysis,
                           cropRasters){

  # get input data
  TimeRange <- as.integer(FHanalysis$TimeSpan)
  TimeNames <- as.character(FHanalysis$TimeSpan)
  LTR <- length(TimeRange)
  r <- FHanalysis$FH_IDr
  FH_ID <- raster::values(r)
  # make raster template
  r <- raster::raster(nrows = nrow(r),
              ncols = ncol(r),
              ext = extent(r),
              crs = crs(r),
              vals = NULL
              )
  # clean memory
  gc()

  # get combination input data
  PLM <- cropRasters$PLM
  EFG <- cropRasters$EFG
  FIRE_REG <- as.integer(cropRasters$RGN)
  FIREFMZ <- as.integer(cropRasters$FIREFMZ)
  DELWP <- cropRasters$DELWP
  #PUBLIC<-as.integer()                                                   #####-----remove
  # combine into data table
  AllCombs<-data.table::as.data.table(cbind(FH_ID,EFG,FIRE_REG,FIREFMZ,PLM,DELWP))
  #TFI<-left_join(EFG_DF,TFI_LUT)                                         #####-----remove
  # remove variables no longer needed
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


  # get the number of pixels in each combination(?)                       #####-----is 'each combination' right?
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
  U_AllCombs_TFI$Hectares <- U_AllCombs_TFI$nPixel * cellsToHectares()

  # return function data
  return(list("U_AllCombs_TFI" = U_AllCombs_TFI,"Index_AllCombs" = Index_AllCombs))
}

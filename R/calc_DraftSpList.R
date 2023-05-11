#' Calculate a draft species list for defined polygon
#' @details Calculate the proportion of cells for the HDM in the region for each species
# works by using the indices of the standard dimensions of raster that are in
# the supplied shapefile region boundary.  The ouptut table contains all species for which there are HDMs. It
#'is intended only as a starting point and requires manual quality control to
#'produce a useful species list for the area by editing the resulting .csv
#'file
#' @param REG_NO integer DELWP fire region number 1:6 ,99 for Statewide analysis, or 7 for ad-hoc boundary polygon default =7 (see look up table REG_LUT for values)
#' @param RasterRes integer 225 - raster resolution is always 225 for this function for speed
#' @param PUBLIC_LAND_ONLY logical whether to restrict analysis to public land only or the whole polygon
#' @param myPoly default clipPoly sf polygon data frame of LF_REGIONs (default) or ad hoc polygon - used in conjunction with REG_NO
#' @param generalRasterDir relative path to directory containing rasters of FIRE_REG, and PUBLIC LAND (PLM_GEN)
#' @param TaxonList path to default species attribute table coantianing paths to HDMs in column HDMPath
#'
#' @return data.frame created from  TaxonList with columns appended for:
#' \itemize{
#' \item cellsInState count of the number of cells in the state within the Binary HDM for the species
#' \item cellsInArea count of the number of cells within myPoly and within the Binary HDM for the species
#' \item areaProp proportion of binary HDM for the state within myPoly
#' }
#' @export
calc_DraftSpList <- function(REG_NO,
                            RasterRes= 225,
                            PUBLIC_LAND_ONLY,
                            myPoly = clipPoly,
                            generalRasterDir = "./InputGeneralRasters",
                            TaxonList = "./ReferenceTables/DraftTaxonListStatewidev2.csv"
                            #,myHDMVals = "./HDMS/HDMVals225.qs"
                            ){
  #HDMVals<-qs::qread(myHDMVals)#loads matrix of binary thrsholded  HDMVals  called "HDMvals"
  REG_NO <- as.integer(as.numeric(REG_NO))                                  ######-----go straight to int? rather than wrap the numeric -did not seem to work when tried as.integer(REG_NO)
  TaxonList <- readr::read_csv(TaxonList)
  CropDetails <- cropToOutput (REG_NO = REG_NO,
                               myRasterRes = RasterRes,
                               PUBLIC_LAND_ONLY = PUBLIC_LAND_ONLY,
                               myPoly = myPoly,
                               generalRasterDir = "./InputGeneralRasters"
  )

  for (i in 1:nrow(TaxonList)){
    if (file.exists(TaxonList$HDMPath[i])){
    r <- terra::rast(TaxonList$HDMPath[i])
    TaxonList$cellsInArea[i] <- sum(terra::extract(x = r,y = CropDetails$inCells),na.rm=T)[1]
    TaxonList$cellsInState[i] <- sum(terra::extract(x = r,y = terra::cells(r)))
    }


  }
  # calc proportion of statwide population
  TaxonList <-TaxonList %>% dplyr::mutate(areaProp = signif(cellsInArea / cellsInState, digits = 2))

  return(TaxonList)
}

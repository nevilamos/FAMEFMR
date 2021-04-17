#'Calculate the species in each EFG in given area for GSO calculations.
#'
#'works by using the indices of the standard dimensions raster that are in the
#'supplied shapefile region boundary (via function cropNAborder )

#' @param REG_NO integer DELWP fire region number 1:6 ,99 for Statewide analysis, or 7 for ad hoc boundary polygon default =7 (see look up table REG_LUT for values)
#' @param RasterRes integer 225 - raster resolution is always 225 for this function for speed
#' @param PUBLIC_LAND_ONLY logical whether to restrict analysis to public land only or the whole polygon
#' @param myPoly default clipPoly sf polygon data frame of LF_REGIONs (default) or ad hoc polygon - used in conjunction with REG_NO
#' @param generalRasterDir relative path to directory containing rasters of FIRE_REG, and PUBLIC LAND (PLM_GEN)
#' @param splist path to default species attribute table default is "./ReferenceTables/DraftTaxonListStatewidev2.csv"
#' @param myHDMVals sparse matrix of cell values for Habitat Distribution Model rasters at 225m pixel size #' saved as a qs file on disk
#' @param TFI_LUT data.frame lookup table for EFG loaded in setup
#' @param myResultsDir path of directory where output will be saved
#' @return list of three data frames LMU_EFG_AREA, Spp_EFG_LMU, and LMU_Scenario used as draft inputs to aspatial GSO calculations
#'@export

calc_Spp_EFG_LMU <- function(REG_NO,
                             RasterRes = 225,
                             PUBLIC_LAND_ONLY,
                             myPoly = clipPoly,
                             generalRasterDir = "./InputGeneralRasters",
                             splist = "./ReferenceTables/DraftTaxonListStatewidev2.csv",
                             myHDMVals = "./HDMS/HDMVals225.qs",
                             myResultsDir= ResultsDir,
                             #EFGRas = EFGRas,
                             TFI_LUT = TFI_LUT){
  # load HDM data
  HDMVals<-qs::qread(myHDMVals)
  mySpList <- utils::read.csv(splist)[,c( "TAXON_ID","COMMON_NAME","SCIENTIFIC_NAME")]
  #get path to correct resolution EFG raster
  EFGRas<-file.path(generalRasterDir,paste0("EFG_NUM_",RasterRes,".tif"))
  #EFG <- raster::values(raster(EFGRas))
  REG_NO <- as.integer(as.numeric(REG_NO))
  CropDetails <- cropNAborder (REG_NO = REG_NO,
                               myRasterRes = RasterRes,
                               PUBLIC_LAND_ONLY = PUBLIC_LAND_ONLY,
                               myPoly = myPoly,
                               generalRasterDir = generalRasterDir
  )

  TFI_LUT<-dplyr::rename(TFI_LUT,EFG_NO = EFG)
  # crop EFG and HDMVals
  EFG <- raster::values(raster::raster(EFGRas))[CropDetails$clipIDX]
  EFG[is.na(EFG)] <- 99
  HDMVals <- HDMVals[CropDetails$clipIDX,]
  mode(EFG) <- "integer"

  # write spp EFG LMU csv
  A <- Matrix.utils::aggregate.Matrix	(HDMVals, EFG ,fun = 'sum')
  myDf <- as.data.frame(as.matrix(A))
  myDf$EFG_NO <- as.integer(rownames(myDf))
  myDf <- tidyr::gather(myDf, key = "TAXON_ID", "CellCount", -EFG_NO)
  myDf <- myDf[myDf$CellCount > 0,]
  myDf$TAXON_ID <- as.integer(myDf$TAXON_ID)
  myDf$EFG_NO <- as.integer(myDf$EFG)
  myDf$ha <- myDf$CellCount * cellsToHectares(RasterMetres = RasterRes)

  myDf <- dplyr::left_join(myDf, TFI_LUT[,c("EFG_NO","EFG_NAME")], by = "EFG_NO")
  myDf <- dplyr::left_join(myDf, mySpList)
  myDf <- myDf[,c("COMMON_NAME","EFG_NO","EFG_NAME","TAXON_ID","CellCount","ha")]

  # write EFG areas csv
  EFG_AREAS <- as.data.frame(table(EFG))
  EFG_AREAS$ha <- EFG_AREAS$Freq* cellsToHectares(RasterMetres = RasterRes)
  EFG_AREAS$EFG_NO <- as.numeric(levels(EFG_AREAS$EFG))
  EFG_AREAS <- dplyr::right_join(TFI_LUT[,c("EFG_NO", "EFG_NAME")], EFG_AREAS)
  EFG_AREAS<-EFG_AREAS[,c("EFG_NO",	"EFG_NAME",	"ha")]
  names(EFG_AREAS)[3]<-"Area"

  #Make template for user to edit for scenarios
  SCENARIO_TEMPLATE <-merge(data.frame("GS_NAME"=c("Juvenile","Adolescent","Mature","Old"),
                                       "GS_ID"=1:4),EFG_AREAS[,c("EFG_NO",	"EFG_NAME")])[,c("EFG_NO",	"EFG_NAME","GS_NAME","GS_ID")]
  SCENARIO_TEMPLATE$Scenario <- "Scenario_0"
  SCENARIO_TEMPLATE$PercLandscape <- 0.25

  return(list ("LMU_EFG_AREA" = EFG_AREAS,"Spp_EFG_LMU" = myDf,"LMU_Scenario"=SCENARIO_TEMPLATE))
}

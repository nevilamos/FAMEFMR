#'Calculate the species in each EFG in given area for GSO calculations.
#'
#'works by using the indices of the standard dimensions raster that are in the
#'supplied shapefile region boundary (via function cropNAborder )

#' @param REG_NO integer DELWP fire region number 1:6 ,99 for Statewide analysis, or 7 for ad hoc boundary polygon default =7 (see look up table REG_LUT for values)
#' @param RasterRes integer 225 - raster resolution is always 225 for this function for speed
#' @param PUBLIC_LAND_ONLY logical whether to restrict analysis to public land only or the whole polygon
#' @param myPoly default clipPoly sf polygon data frame of LF_REGIONs (default) or ad hoc polygon - used in conjunction with REG_NO
#' @param generalRasterDir relative path to directory containing rasters of FIRE_REG, and PUBLIC LAND (PLM_GEN)
#' @param TaxonList path to  species attribute table containing paths to HDMs default is "./ReferenceTables/FAME_TAXON_LIST.csv"
#' @param TFI_LUT data.frame lookup table for EFG loaded in setup
#' @param myResultsDir path of directory where output will be saved
#' @return list of three data frames LMU_EFG_AREA, Spp_EFG_LMU, and LMU_Scenario used as draft inputs to aspatial GSO calculations
#'@export

make_Draft_GSO_inputs <- function(REG_NO,
                             RasterRes = 225,
                             PUBLIC_LAND_ONLY,
                             myPoly = clipPoly,
                             generalRasterDir = "./InputGeneralRasters",
                             TaxonList= "D:/FAMEshiny/ReferenceTables/FAME_TAXON_LIST.csv",
                             myResultsDir= ResultsDir,
                             TFI_LUT = TFI_LUT){
  # load HDM data
  #HDMVals<-qs::qread(myHDMVals)
  mySpList <- readr::read_csv(TaxonList)
  mySpList<-mySpList %>% dplyr::filter(Include == "Yes")
  myHDMs<-mySpList$HDMPath
  #get path to correct resolution EFG raster
  EFGRas<-file.path(generalRasterDir,paste0("EFG_NUM_",RasterRes,".tif"))
  #EFG <- terra::values(raster(EFGRas))
  REG_NO <- as.integer(as.numeric(REG_NO))
  CropDetails <- cropToOutput (REG_NO = REG_NO,
                               myRasterRes = RasterRes,
                               PUBLIC_LAND_ONLY = PUBLIC_LAND_ONLY,
                               myPoly = myPoly,
                               generalRasterDir = generalRasterDir

  )
  if ("EFG" %in% names(TFI_LUT)){

  TFI_LUT<-dplyr::rename(TFI_LUT,EFG_NO = EFG)
  }
  # crop EFG and HDMVals
  myHDMRasters<-terra::rast(myHDMs)
  EFG <- terra::extract(terra::rast(EFGRas),CropDetails$inCells)
  EFG[is.na(EFG)] <- 99
  EFG<-as.matrix(EFG)
  HDMVals <- terra::extract(myHDMRasters,CropDetails$inCells)
  HDMVals<-as.matrix(HDMVals)
  HDMVals[HDMVals>0]<-1
  HDMVals[is.na(HDMVals)]<-0
  dimnames(HDMVals)[[2]]<-mySpList$TAXON_ID
  mode(EFG) <- "integer"

  # write spp EFG LMU csv
  #A <- Matrix.utils::aggregate.Matrix	(HDMVals, EFG ,fun = 'sum')
  aggregate(HDMVals,list(EFG),FUN='sum')
  myDf <- aggregate(HDMVals,list(EFG),FUN='sum')
  names(myDF)[1]<-EFG_NO
  myDf <- tidyr::gather(myDf, key = "TAXON_ID", "CellCount", -EFG_NO)
  myDf <- myDf[myDf$CellCount > 0,]
  myDf$TAXON_ID <- as.integer(myDf$TAXON_ID)
  myDf$EFG_NO <- as.integer(myDf$EFG)
  myDf$ha <- myDf$CellCount * cellsToHectares(RasterMetres = RasterRes)

  myDf <- dplyr::left_join(myDf, TFI_LUT[,c("EFG_NO","EFG_NAME")], by = "EFG_NO")
  myDf <- dplyr::left_join(myDf, mySpList)%>%
    tidyr::drop_na()
  myDf <- myDf[,c("COMMON_NAME","EFG_NO","EFG_NAME","TAXON_ID","CellCount","ha")]

  # write EFG areas csv
  EFG_AREAS <- as.data.frame(table(EFG))
  EFG_AREAS$ha <- EFG_AREAS$Freq* cellsToHectares(RasterMetres = RasterRes)
  EFG_AREAS$EFG_NO <- as.numeric(levels(EFG_AREAS$EFG))
  EFG_AREAS <- dplyr::right_join(TFI_LUT[,c("EFG_NO", "EFG_NAME")], EFG_AREAS)%>%
    tidyr::drop_na()
  EFG_AREAS<-EFG_AREAS[,c("EFG_NO",	"EFG_NAME",	"ha")]
  names(EFG_AREAS)[3]<-"Area"

  #Make template for user to edit for scenarios
  SCENARIO_TEMPLATE <-merge(data.frame("GS_NAME"=c("Juvenile","Adolescent","Mature","Old"),
                                       "GS_ID"=1:4),EFG_AREAS[,c("EFG_NO",	"EFG_NAME")])[,c("EFG_NO",	"EFG_NAME","GS_NAME","GS_ID")]%>%
    tidyr::drop_na()
  SCENARIO_TEMPLATE$Scenario <- "Scenario_0"
  SCENARIO_TEMPLATE$PercLandscape <- 0.25

  return(list ("LMU_EFG_AREA" = EFG_AREAS,"Spp_EFG_LMU" = myDf,"LMU_Scenario"=SCENARIO_TEMPLATE))
}

############################################################################
#Utility functions used in FAME Fire History analyses (often called within the
#main analysis functions) some are likely to be project specific
#others may be useful more generically written by nevil.amos@delwp.vic.gov.au
############################################################################


#' remove empty directories from path
#' @details Removes all empty subdirectories from the nominated path does not remove the nominated path directory even if empty
#' @param rootDir relative path to remove all empty subdirectories from. Default
#'   value "./Results"
#' @export
removeEmptyDirs <- function(rootDir="./Results"){

  dirList <- list.dirs(rootDir)[-1] #[-1] makes sure the root results directory is not deleted
  if(length(dirList) > 0){
    delList <- dirList[sapply(dirList, function(x) length(list.files(x)) == 0)]
    while((length(delList) > 0) & (length(dirList) > 0)){
      unlink(delList, recursive = TRUE)
      dirList <- list.dirs(rootDir)[-1]
      delList <- dirList[sapply(dirList, function(x) length(list.files(x)) == 0)]
    }
  }
}


#' Adds concatenated String of X and Y coordinates of centroids of polygons to
#' Simple Features polygon object. This String acts as a key to identify
#' spatially identical polygons for use in tidyverse pivot functions.
#' @param myDF sf polygon object
#' @return character vector of XYStrings
#' @export
add_xystring <- function(myDF){
  Coords <- as.data.frame(
    do.call(
      rbind,
      sf::st_geometry(sf::st_centroid(myDF))))
  names(Coords) <- c("X","Y")
  XYString <- paste(Coords$X, Coords$Y, sep = "")
  x <- dplyr::mutate(myDF, XYString)
  return(x)
}

#' Extract VBA (Victorian Biodiversity Atlas)  species ID numbers from file
#' paths extracts four or five digit species numbers (Victorian Biodiversity
#' Atlas TAXON_IDs) from vector of paths or file names containing files of e.g.
#' species HDMS containing the 5 digit TAXON_ID in their name
#' @param x Vector of species file Pathnames containing VBA numbers
#' @return numeric vector of 4or 5 digits (ususally TAXON_ID)
#' @export
get_Spp_No <- function( x = "Vector of Sp file Pathnames"){
  Fnames= basename(x)
  pos = regexpr("[0-9][0-9][0-9][0-9]+", Fnames)
  myEnd = pos-1+attr(pos, "match.length")
  y = as.numeric(substr(Fnames, pos, myEnd))
  return(y)
}

#' Joins one or more lookup tables to table containing ID values Function joins
#' Lookup tables (LUTS) to dataframe containing ID_NO: Name combinations
#'
#' @param myDF dataframe or similar containing indices for the LUTS listed, to
#'   which the LUTS will be dplyr::left_joined
#' @param LUTS vector of names of LUTS in memory defaults
#'   =c("TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")
#' @return a data.frame with the LUTS joined to it
#' @export
Join_Names <- function(myDF,   #dataframe or similar containing indices for the LUTS listed
                       LUTS = c("TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")){
  for(i in LUTS){
    try(myDF <- dplyr::left_join(myDF,get(i)))}
  return (myDF)
}


#' Calculates multiplier to convert from raster cell count to area in hectares
#' @param RasterMetres numeric Value cell resolution in Metres (usually from
#'   RasterRes in settings file).
#' @return numeric Multiplier to convert cell count to area in hectares
#' @export
cellsToHectares <- function(RasterMetres = RasterRes){
  (RasterMetres / 100) ^ 2
}


#' Checks whether all values in  one vector are in another vector
#' @param x Vector of values to check if all are in second vector
#' @param v Second vector of values that may or may not contain all values in x
#' @return logical
#' @export
notAllIn <- function(x, v = V){
  anyNA(unlist(lapply(x, match, v)))
}



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

#' Fix Pivot_wider list of lists columns
#' @details Supporting function to deal with pivot_wider returning list of lists in some cases
#' @param df wide format data frame with fields that are lists of lists
#'
#' @return wide format data frame without fields that are lists of lists
#' @export

unlistPivot_wider <- function(df){
  df[unlist(lapply(df , is.null))] <- NA
  Y <- unlist(df)
  return (Y)
}



#' Makes Long format  Fauna  abundance table
#'
#' @details Supporting function to Make long format table for fauna abundance scores
#'  by "VBA_CODE"  ,FireType EFG and Growth Stage from input wide format table currently
#'   deals only with FireType of "High" and "Low"
#'   which are converted to 2 and 1 respectively
#' @param AbundDataByGSFile .csv input file containing fields
#' "EFG_NO", "GS4_NO",  "FireType" , "Abund", "VBA_CODE"
#' with Abund values for' both FireTypes for each growth stage "GS4_NO"
#' @param myEFG_TSF_4GS table of each combination of "EFG_NO", "GS4_NO", and "YSF" generally read in at beginning of FAME in settings file
#' @return long format table with one row for each combination of "EFG_NO", "GS4_NO",  "FireType" , "Abund", "VBA_CODE" and "YSF"
#' sorted by VBA_CODE
#' @export

makeAbundDataLong<- function(AbundDataByGSFile = "./ReferenceTables/OrdinalExpertLong.csv",
                             myEFG_TSF_4GS = EFG_TSF_4GS){
  AbundDataByGS <- utils::read.csv(AbundDataByGSFile)[,c("EFG_NO", "GS4_NO",  "FireType" , "Abund", "VBA_CODE")]
  AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="High"]<-2
  AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="Low"]<-1
  AbundDataByGS<-AbundDataByGS[!is.na(AbundDataByGS$Abund),c("EFG_NO", "GS4_NO",  "FireTypeNo" , "Abund", "VBA_CODE")]


  AbundDataLong = merge(AbundDataByGS, myEFG_TSF_4GS,   by=c('EFG_NO','GS4_NO'))
  AbundDataLong<-AbundDataLong[order(AbundDataLong$VBA_CODE),]
  return(AbundDataLong)
}





## FUNCTION makeHDMValsfromRasters  ------------------------------------------------------
#' Function makes a matrix of HDM values(1,NA) constrained to those cells that
#' are indexed in the cropped area
#' @param myHDMSpp_NO vector of TAXON_IDs for species to be included in output
#' @param myCropRasters list of rasters and indices and cell values created by function cropNAborder()
#' @importFrom iterators iter
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @export
makeHDMValsfromRasters <- function(myHDMSpp_NO = HDMSpp_NO,
                                   myCropRasters = cropRasters)
{
  HDMPaths <- dir(myCropRasters$HDM_RASTER_PATH, full.names=TRUE, pattern = ".tif$")
  HDMPaths <- HDMPaths[get_Spp_No(HDMPaths) %in%
                         myHDMSpp_NO]

  print("reading HDMvalues")

  # set up parallel processing
  cl<-parallel::makeCluster(Ncores)
  doParallel::registerDoParallel(cl, cores=Ncores)

  myHDMVals <- foreach::foreach(i = iterators::iter(HDMPaths),.combine = cbind,.packages = "raster") %dopar% {
    myVals <- raster::values(raster::raster(i))[myCropRasters$IDX]
    myVals}
  parallel::stopCluster(cl)

  colnames(myHDMVals) <- as.character(get_Spp_No(HDMPaths))
  return(myHDMVals)
}




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



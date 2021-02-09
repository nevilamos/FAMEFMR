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
#'   which the LUTS will be left_joined
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

#' Calculate a draft species list for defined polygon
#' @details Calculate the proportion of cells for the HDM in the region for each species
# works by using the indices of the standard dimensions of raster that are in
# the supplied shapefile region boundary.  The ouptut table contains all species for which there are HDMs. It
#'is intended only as a starting point and requires manual quality control to
#'produce a useful species list for the area by editing the resulting .csv
#'file
#' @param REG_NO integer DELWP fire region number 1:6 ,99 for Statewide analysis, or 7 for ad hoc boundary polygon default =7 (see look up table REG_LUT for values)
#' @param RasterRes integer 225 - raster resolution is always 225 for this fucntion for speed
#' @param PUBLIC_LAND_ONLY logical whether to restrict analysis to public land only or the whole polygon
#' @param myPoly default clipPoly sf polygon data frame of LF_REGIONs (default) or ad hoc polygon - used in conjunction with REG_NO
#' @param generalRasterDir relative path to directory containing rasters of FIRE_REG, and PUBLIC LAND (PLM_GEN)
#' @param splist path to default species attribute table default is
#' @param HDMVals
#'
#' @return data.frame created from splist with columns appended for:
#' \itemize{
#' \item cellsInState count of the number of cells in the state within the Binary HDM for the species
#' \item cellsInArea count of the number of cells within myPoly and within the Binary HDM for the species
#' \item areaProp proportion of binary HDM for the state within myPoly
#' }
#' @export
calcDraftSpList <- function(REG_NO,  # can this match look up table REG_LUT     #REG_NO of defined region from input (1:6) or 99 for statewide or 7 for Ad Hoc Poly),
                            RasterRes= 225, # raster resolution of 225is used for this function for speed - it is only required to get approximate values so finer resolution in not needed
                            PUBLIC_LAND_ONLY,
                            myPoly = clipPoly,       #shapefile of LF_REGIONs (default) or adhoc region,   ######-----should this be clipPoly for consistency?
                            generalRasterDir = "./InputGeneralRasters",
                            splist = "./ReferenceTables/DraftTaxonListStatewidev2.csv",#always uses the default species list
                            HDMVals = HDMVals225){                          ######-----HDMVals225 is commented out in settings, and not called anywhere prior.
  load(HDMVals)
  REG_NO <- as.integer(as.numeric(REG_NO))                                  ######-----go straight to int? rather than wrap the numeric -did not seem to work when tried as.integer(REG_NO)
  splist <- read.csv(splist)
  CropDetails <- cropNAborder (REG_NO = REG_NO,
                               RasterRes = RasterRes,
                               PUBLIC_LAND_ONLY,
                               myPoly = myPoly,
                               generalRasterDir = "./InputGeneralRasters"
  )

  # get cells in polygon
  cellsInArea <- colSums(HDMVals[CropDetails$clipIDX,])
  cellsInState <- colSums(HDMVals)

  # calc proportion of statwide population
  areaProp <- signif(cellsInArea / cellsInState, digits = 2)
  # pull data into dataframe
  TAXON_ID <- as.numeric(colnames(HDMVals))
  myDF <- data.frame(TAXON_ID, cellsInState, cellsInArea, areaProp)
  myDF <- left_join(splist, myDF)
  return(myDF)
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

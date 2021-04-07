############################################################################
#Utility functions used in FAME Fire History analyses (often called within the
#main analysis functions) written by nevil.amos@delwp.vic.gov.au
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


#' Calculate Geometric mean
#' @details Calculates geometric mean of a vector of positive numeric values

#' @param x numeric vector
#'
#'
#' @return numeric or NA if values are not all > 0
geoMean<-function (x){
  if(any(is.na(x))){
    y = NA
    warning("vector contains NA values")
    } else if(any(x<=0)){
      y=NA
      warning("Not all x are positive values")
    } else {
        y=prod(x)^(1/length(x))
        }
  return(y)
}













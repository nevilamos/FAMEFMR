############################################################################
#Utility functions used in FAME Fire History analyses ( often called within the
#main analysis functions) some are likely to be only of project specific use
#others may be useful more generically written by nevil.amos@delwp.vic.gov.au
############################################################################


## FUNCTION TO REMOVE ALL EMPTY DIRECTORIES IN A PATH ----------------------
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

## ADD XY CENTROID STRING --------------------------------------------------
# Adds single String of X and Y of centroids of polygons to Simple Features polygon object
# This String acts as a key to identify spatially identical polygons for use in tidyverse pivot functions

add_xystring <- function(myDF){
  require(dplyr)
  Coords <- as.data.frame(do.call(rbind,st_geometry(st_centroid(myDF))))
  names(Coords) <- c("X","Y")
  XYString <- paste(Coords$X, Coords$Y, sep = "")
  x <- mutate(myDF, XYString)
  return(x)
}

## GET VIC SPP ID No. FROM HDM PATHS ---------------------------------------
#extracts four or five digit species numbers (Victorian Biodiversity Atlas
#TAXON_IDs) from  basenames of paths containing files of e.g. species HDMS
#containing the 5 digit TAXON_ID in their name
get_Spp_No <- function(x = "Vector of Sp file Pathnames"){
  Fnames= basename(x)
  pos = regexpr("[0-9][0-9][0-9][0-9]+", Fnames)
  myEnd = pos-1+attr(pos, "match.length")
  y = as.numeric(substr(Fnames, pos, myEnd))
  return(y)
}

## FUNCTION Join_Names ------------------------------------------------------
# Function joins Lookup tables to dataframe containing ID_NO: Name combinations
# at the moment the LUTS are hard wired.

Join_Names <- function(myDF,   #dataframe or similar containing indices for the LUTS listed 
                       LUTS = c("TFI_LUT","TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")){
  for(i in LUTS){
    try(myDF <- left_join(myDF,get(i)))}
  return (myDF)
}


## FUNCTION cellsToHectares -------------------------------------------------
# Function to calclulate multiplier from cells to hecatres
# cell resolution in Metres (usually from RasterRes in settings file)
cellsToHectares <- function(RasterMetres = RasterRes){
  (RasterMetres / 100) ^ 2
}


## FUNCTION notAllIN -------------------------------------------------------
# function to check whether all values in vector are in another vector of
# permitted values, returns  TRUE/FALSE
notAllIn <- function(x, v = V){
  anyNA(unlist(lapply(x, match, v)))
}
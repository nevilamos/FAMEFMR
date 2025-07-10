#' Wrapper for to complete fhProcess1 and fhProcess2
#'
#' @details The function takes a shapefile or  geopackage layer of Fire history
#'   containing polygons with two fields: FIRETYPE and SEASON Where polygons of
#'   different FIRETYPE or SEASON overlap the function constructs unique
#'   non-overlapping polygon of their intersections ( and non intersecting areas
#'   ) and attributes each polygon with sequential fire SEASON (SEAS01, SEAS02
#'   ...) and corresponding FIRETYPE (TYPE01,TYPE02 ...)
#' @details It then calculates all the intervals between sequential fires, and
#'   Time Since fire (TSF) and Last Fire Type (LFT) and Last burnt year (LBY)
#'   for each SEASON as defined in the input arguments, these values are append
#'   to the output sf polygon dataframe.
#' @param firstFH path to the input fire history geodatabase (".gdb" or ".gpkg")
#'   or shapefile (".shp") usually  provided in settings
#' @param firstFHLayer if rawFH is a geodatabase the name of the layer containing
#'   the Fire History if this is not provided and the raw FH is a geodatabase or
#'   geopackage then the first layer in the file is returned.
#' @param secondFH Second fire history to be combined with FH1 to make a fire
#'  scenario same formats as for firstFH or NULL if there is only a single FH used.
#' @param secondFHLayer if secondFH is a geodatabase the name of the layer containing
#'   the Fire History if this is not provided and the raw FH is a geodatabase or
#'   geopackage then the first layer in the file is returned.
#' @param start.SEASON integer First SEASON for which output is wanted (four
#'   digit year as integer), if NUll then second season in in history is used
#'   (cannot use first season because it has no interval, this may still fail if
#'   there is no overlap)
#'
#' @param end.SEASON  integer Last SEASON required, if NULL then largest value
#'   in fire history scenario used
#'
#' @param OtherAndUnknown integer Value to use for cases where fire type is:
#'   "OTHER" or "UNKNOWN" = NA, "BURN" = 1, "BUSHFIRE" = 2. NA = Fire excluded
#'   from analysis. usually set in settings file
#'
#' @param validFIRETYPE vector of valid names in the input FIRETYPE column in
#'   the input fire history dataset(s), if the column contains NA or values not
#'   on this list an error will occur
#' @param max_interval = NULL #integer number of years  maximum inter fire
#'   interval after a HIGH fire for which subsequent fire is reassigned to HIGH
#'   default is NULL in which case there is no reassignment
#' @param baseFire Default NULL otherwise four digit integer SEASON for fire
#'   applied across the whole bounding box
#'
#' @return A list containing: \itemize{
#' \item OutDF sf polygons dataframe containing all the fire history attributes
#' \item TimeSpan integer vector sequence of SEASONS to in the analysis output
#' \item YSFNames names of TSF years in output, needed by downstream functions
#' \item LBYNames names of LBY years in output, needed by downstream functions
#' \item LFTNames names of LBY years in output, needed by downstream functions }
#'
#' @export
#'

fhProcess<-function(firstFH,
                    firstFHLayer=NULL,
                    secondFH=NULL,
                    secondFHLayer=NULL,
                    OtherAndUnknown=2,
                    baseFire=baseFire,
                    max_interval = 0,
                    start.SEASON = NULL,
                    end.SEASON = NULL,
                    validFIRETYPE = c("BURN","BUSHFIRE","OTHER","UNKNOWN")
){


  outFH1<-fhProcess1(inFH = firstFH,inFHLayer=firstFHLayer,secondFH=secondFH,secondFHLayer=secondFHLayer,OtherAndUnknown=OtherAndUnknown,baseFire=baseFire,validFIRETYPE = validFIRETYPE)
  outFH2<-fhProcess2(inFH1 = outFH1,
                              start.SEASON=start.SEASON,
                              end.SEASON=end.SEASON,
                              max_interval =max_interval
  )
  return(outFH2)
}



#' Main Fire History Fire Sequence analysis function
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
#'
#' @param rawFH path to the input fire history geodatabase (".gdb" or ".gpkg")
#'   or shapefile (".shp") usually  provided in settings
#' @param rawFHLayer if rawFH is a geodatabase the name of the layer containing
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
#' @param validFIRETYPE character vector of valid FIRETYPE values for checking
#'   the input file , provided in settings file.
#'
#' @return A list containing: \itemize{
#' \item OutDF sf polygons dataframe containing all the fire history attributes
#' \item TimeSpan integer vector sequence of SEASONS to in the analysis output
#' \item YSFNames names of TSF years in output, needed by downstream functions
#' \item LBYNames names of LBY years in output, needed by downstream functions
#' \item LFTNames names of LBY years in output, needed by downstream functions }
#' @export
fhProcess<-function(rawFH = "path of the rawFH file geopackage or gdb to use",
                    rawFHLayer = NULL,
                     start.SEASON = NULL,    # first season for which output is wanted (four digit year as integer), if NUll then second season in in history is used (cannot use first season because it has no interval, this may still fail if there is no overlap)
                     end.SEASON = NULL,      # last season required, if NULL then largest value in fire history scenario used
                     OtherAndUnknown,     # ## link to look up table FIRETYPE_LUT?? ## Default is 2 ("BUSHFIRE").  (2,1,NA) value to use for cases where fire type is: "OTHER" or "UNKNOWN" = NA, "BURN" = 1, "BUSHFIRE" = 2. NA = Fire excluded from analysis.
                     #####----- the OtherAndUnknown default should be NA? currently you include bushfires (unless otherwise stated)
                     validFIRETYPE =c("BURN", "BUSHFIRE", "UNKNOWN", "OTHER")
){
  # read in shapefile

  # error checks------
  if(tools::file_ext(rawFH)%in%c("shp","gdb","gpkg")){
  if( is.null(rawFHLayer)){
    myDF <- sf::st_read(dsn = rawFH)
  }else{
    myDF <- sf::st_read(dsn = rawFH,layer = rawFHlayer)
  }
  }else{stop("rawFH file is not a shapefile,geopackage or ESRI geodatabase")}




  # check that the input fire History (rawFH) contains the two required fields
  # and that these do not have a missing values

  myDFNames <- names(myDF)
  if(!"sf"%in%class(myDF)){stop("rawFH file is not a spatial dataset")}
  if (!st_geometry_type(myDF,by_geometry = F)[1]%in%c("POLYGON","MULTIPOLYGON")){stop("rawFH file is not a polygon dataset")}
  if((!"SEASON" %in% myDFNames)) stop ("rawFH file does not contain field 'SEASON'")
  if(!class(myDF[["SEASON"]]) %in% c("numeric","integer"))stop ("rawFH file field 'SEASON' is not type integer")
  if(anyNA(myDF$SEASON)) stop ("rawFH file has missing values in  field 'SEASON'")
  if((!"FIRETYPE"%in% myDFNames)) stop ("rawFH file does not contain field 'FIRETYPE'")
  if(notAllIn(x = myDF$FIRETYPE, v = validFIRETYPE)) stop ("rawFH shapefile has missing or invalid values in field 'FIRETYPE'")



  # timespan (range of consecutive years) for which Fire History sequences are calculated
  # from start and end seasons defined in the settings file
  # else if NULL; start and end season for calculation is determined by the min+1 and max values in the dataset
  # the earliest start season is the date of a second season in the input data (min+1), next 3 lines prevent manual start
  # season setting less than this value, which would cause an error.
  min.SEASON<-sort(unique(myDF$SEASON))[2] # second season in fire history
  if(is.null(start.SEASON)){
    start.SEASON = min.SEASON
  } else {
    if(start.SEASON < min.SEASON){
      start.SEASON = min.SEASON
    } else {
      start.SEASON = start.SEASON
    }
  }

  if(is.null(end.SEASON)){
    max.SEASON <- max(myDF$SEASON)
  } else {
    max.SEASON = end.SEASON
  }

  TimeSpan<-start.SEASON:max.SEASON

  # recode FIRETYPE to integer values value for Other and Unknown is set in settings file   #####-----can the look up table be used?
  myDF$FIRETYPE_NO[myDF$FIRETYPE == "BURN"] <- 1
  myDF$FIRETYPE_NO[myDF$FIRETYPE == "BUSHFIRE"] <- 2
  myDF$FIRETYPE_NO[myDF$FIRETYPE == "OTHER"] <- OtherAndUnknown
  myDF$FIRETYPE_NO[myDF$FIRETYPE == "UNKNOWN"] <- OtherAndUnknown            #####-----currently set to 2 (bushfire)


  #this line was added to force myDF$FIRETYPE_NO to integer value - otherwise in shiny version of function was becoming character
  myDF$FIRETYPE_NO<-as.integer(myDF$FIRETYPE_NO)

  # add a field containing a string representation of the x and y values of the
  # centroid of each polygon in myDF to use as a unique identifier for unique
  # polygons in the FH analysis (ie two spatially identical polygons with
  # different SEASON of FireType values will share an unique xycoords string) this
  # is used to get the wide format sequence for each spatially unique polygon
  myDF <- add_xystring(myDF)


  # filter spatial dataframe by unique combinations polygons geometry and XYString
  #with SEASON and FIRETYPE

  myDF<-unique(myDF[,names(myDF)]
  )

  # order the spatially unique polygons by SEASON then FIRETYPE_NO
  # adds a sequential number to them, this is what allows the flattening and subsequent
  # reduction to fire sequences without an intervening no fire year
  myDF<-myDF[with(myDF,order(XYString, SEASON, FIRETYPE_NO)),]

  myDF$Sequence<-1
  myDF$Sequence<-unlist(lapply(split(myDF$Sequence, myDF$XYString), cumsum))

  # reduction to the unique sequences of fires with no gaps using dplyr::pivot_wider
  # in order to acheive this efficently the sf object made from the rawFH file is converted
  # to a standard data frame with a geometry column then converted back to an sf object at the end of the process
  print("Making OutDF")
  OutDF <- myDF %>%
    dplyr::select(XYString, Sequence, FireType = FIRETYPE_NO, SEAS = SEASON) %>%
    dplyr::mutate(Sequence = sprintf("%02d", Sequence)) %>%
    as.data.frame %>%
    tidyr::pivot_wider(names_from = Sequence, values_from = c(FireType, SEAS), names_sep="")

  # get just the names for the season and FireType sequence columns
  SEASNames <- names(OutDF)[grep(pattern = "SEAS", names(OutDF))]
  FTNames <- names(OutDF)[grep(pattern = "FireType", names(OutDF))]

  # calculates all inter-fire intervals note by offsetting matrices by one column
  SEAS_Matrix <- as.matrix(OutDF[, SEASNames])
  FT_matrix <- as.matrix(OutDF[, FTNames])
  Cols <- ncol(SEAS_Matrix)
  SEAS_Matrix[SEAS_Matrix == 0] <- NA
  Interval <- SEAS_Matrix[,2:Cols] - SEAS_Matrix[,1:Cols-1]
  IntNames <- paste("INT", sprintf("%02d", 1:(Cols-1)), sep = "")
  colnames(Interval) <- IntNames

  # bind the intervals calcuated back to the main dataframe
  OutDF <- cbind(OutDF, Interval) %>%
    dplyr::select(-XYString)

  OutDF <- sf::st_as_sf(OutDF)


  # calculating the last burnt season for each sequences for each year
  # this process is duplicated here and in CALC_TFI2 function
  # ideally it should only be run once but it is pretty fast so probably does not matter too much.
  LTR <- length(TimeSpan)
  SEAS_Matrix[is.na(SEAS_Matrix)] <- 0
  LBY <- matrix(NA, nrow(SEAS_Matrix), LTR)
  for(i in 1:LTR){
    try({
      y = TimeSpan[i]
      LBY[,i] <- LBY_f(M = SEAS_Matrix, y)
      #print(y)
    })
  }

  # subtraction leads to matrix transpostion which needs to be retransposed in the next couple of lines
  tYSF <- TimeSpan-t(LBY)
  YSF <- t(tYSF)

  # name matrix columns
  YSFNames <- paste0("YSF", TimeSpan)
  LBYNames <- paste0("LBY", TimeSpan)
  LFTNames <- paste0("LFT", TimeSpan)
  colnames(YSF) <- YSFNames
  colnames(LBY) <- LBYNames

  # Create matrix used for getting last firetype by year
  print("calculating lookup matrix for getting last FireType by SEASON")
  SEAS_Matrix[SEAS_Matrix == 0] <- NA
  LUM <- matrix(NA, nrow(SEAS_Matrix), max.SEASON)
  for (i in 1:nrow(SEAS_Matrix)){
    R <- i
    C <- as.numeric(stats::na.omit(SEAS_Matrix[i,]))
    V <- (FT_matrix[i,(1:length(C))])
    LUM[R,C] <- V
  }


  # Calculate last fire type
  print("calculating last fire type")
  LFT <- matrix(NA, nrow(SEAS_Matrix), LTR)
  for(i in 1:nrow(SEAS_Matrix)){
    LFT[i,] <- LUM[i, LBY[i,]]
  }
  colnames(LFT) <- LFTNames


  # put together output dataframe
  OutDF <- cbind(OutDF, YSF)
  OutDF <- cbind(OutDF, LBY)
  OutDF <- cbind(OutDF, LFT)
  OutDF$ID <- as.integer(rownames(OutDF))
  print("Completed making FH object")

  # clean memory
  gc()

  # compile output and return
  results <- list(OutDF = OutDF, TimeSpan = TimeSpan, YSFNames = YSFNames,
                  LBYNames = LBYNames, LFTNames = LFTNames)
  return(results)
}

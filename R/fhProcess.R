#' complete Fire History Processing using QGIS_process
#'
#' @details The function takes a shapefile or  geopackage layer of Fire history
#'   containing polygons with two fields: FIRETYPE and SEASON Where polygons of
#'   different FIRETYPE or SEASON overlap the function constructs unique
#'   non-overlapping polygon of their intersections ( and non intersecting areas
#'   ) and attributes each polygon with sequential fire SEASON (SEAS01, SEAS02
#'   ...) and corresponding FIRETYPE (TYPE01,TYPE02 ...)
#' @details It then calculates all the intervals between sequential fires, and
#'   Time Since fire (TSF) and Last Fire Type (LFT) and Last burnt year (LBY)
#' for each SEASON as defined in the input arguments, these values are append to
#' the output sf polygon dataframe.
#' @param firstFH path to the input fire history geodatabase (".gdb" or ".gpkg")
#'   or shapefile (".shp") usually  provided in settings
#' @param firstFHLayer if rawFH is a geodatabase the name of the layer
#'   containing the Fire History if this is not provided and the raw FH is a
#'   geodatabase or geopackage then the first layer in the file is returned.This
#'   argument is available only in the r script version - it is not exposed in
#'   the shinyApp where it is set to NULL
#' @param secondFH Second fire history to be combined with FH1 to make a fire
#'   scenario same formats as for firstFH or NULL if there is only a single FH
#'   used.
#' @param secondFHLayer if secondFH is a geodatabase the name of the layer
#'   containing the Fire History if this is not provided and the raw FH is a
#'   geodatabase or geopackage then the first layer in the file is returned.This
#'   argument is available only in the r script version - it is not exposed in
#'   the shinyApp where it is set to NULL
#' @param clipShape vector polygon(s)as file path or name of sf vector object in
#'   layer providing the boundary of the area of interest for the analysis. The
#'   fire history inputs and output will be clipped to this polygon layer.
#'
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


fhProcess <- function(firstFH,
                      firstFHLayer = NULL,
                      secondFH = NULL,
                      secondFHLayer = NULL,
                      clipShape = NULL,
                      start.SEASON = NULL,
                      end.SEASON = NULL,
                      OtherAndUnknown = 2,
                      max_interval = 0,
                      validFIRETYPE = c("BURN", "BUSHFIRE", "OTHER", "UNKNOWN"),
                      baseFire = baseFire,
                      precsision =1) {
  print("performing FH processing using qgis native union")
  #check FH inputs and initial combination and cropping of FH files----

  #if more than one FH input combine them

  mySF <-
    fhCheck(inFH = firstFH, inFHLayer = firstFHLayer, validFIRETYPE)
  if (!(is.null(secondFH) | length(secondFH) == 0)) {
    mySF2 <-
      FAMEFMR::fhCheck(inFH = secondFH, inFHLayer = secondFHLayer, validFIRETYPE)
    mySF <- dplyr::bind_rows(mySF, mySF2)
  }



  #calculate bounding box for FH and add base fire for bounding box if required
  BBOX <-
    sf::st_cast(sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(mySF))),
                "POLYGON")
  names(BBOX)[1] <- "geometry"
  if (is.null(baseFire)) {

  }
  else if (nchar(baseFire == 4)) {
    BBOX <-
      sf::st_as_sf(dplyr::bind_cols(BBOX, data.frame(
        SEASON = baseFire,
        FIRETYPE = "BUSHFIRE"
      )))

    mySF <- rbind(BBOX, mySF[, c("SEASON", "FIRETYPE")])

  }
  else {
    stop("baseFire argument is invalid")
  }



  # #crop to area of interest using Terra since faster it also allows this clip to
  # #handle a file path or an sf object name as input
  # if (!is.null(clipShape)) {
  #   mySF <- st_as_sf(terra::crop(vect(mySF), vect(clipShape)))
  # }
  mySF<-result <- qgisprocess::qgis_run_algorithm(
    "native:clip",
    INPUT = mySF,
    OVERLAY = clipShape
    ) %>% sf::st_as_sf()

  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BURN"] <- 1
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BUSHFIRE"] <- 2
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "OTHER"] <- OtherAndUnknown
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "UNKNOWN"] <- OtherAndUnknown



  #perform dissolve and union using qgis_process algorithm ----
  #to produce set of
  #distinct fire history polygons that only overlap with other polygons of the
  #same geometry but different SEASON or FIRETYPE_NO

  FH_diss <- qgisprocess::qgis_run_algorithm(
    "native:dissolve",
    INPUT = mySF,
    FIELD = c("FIRETYPE", "SEASON"),
    SEPARATE_DISJOINT = FALSE
  ) %>%
    st_as_sf()
  FH_union <- qgisprocess::qgis_run_algorithm("native:union",
                                              INPUT = FH_diss,
                                              GRID_SIZE = precsision) %>%
    st_as_sf() %>% filter(st_geometry_type(.) %in%c( "MULTIPOLYGON","POLYGON"))



  #main FHprocess flattening above  -----
  myDF <- add_xystring(FH_union)
  min.SEASON <- sort(unique(myDF$SEASON))[2]
  if (is.na(start.SEASON)) {
    start.SEASON = min.SEASON
  }
  else {
    if (start.SEASON < min.SEASON) {
      start.SEASON = min.SEASON
    }
    else {
      start.SEASON = start.SEASON
    }
  }
  if (is.na(end.SEASON)) {
    max.SEASON <- max(myDF$SEASON)
  }
  else {
    max.SEASON = end.SEASON
  }
  TimeSpan <- start.SEASON:max.SEASON
  myDF <- unique(myDF[, names(myDF)])
  myDF <- myDF[with(myDF, order(XYString, SEASON, FIRETYPE_NO)),]
  myDF$Sequence <- 1
  myDF$Sequence <-
    unlist(lapply(split(myDF$Sequence, myDF$XYString),
                  cumsum))
  print("Making OutDF")
  OutDF <-
    myDF %>% dplyr::select(XYString, Sequence, FireType = FIRETYPE_NO,
                           SEAS = SEASON) %>%
    dplyr::mutate(Sequence = sprintf("%02d",Sequence)) %>%
    as.data.frame %>% tidyr::pivot_wider(names_from = Sequence,
                                         values_from = c(FireType, SEAS),
                                         names_sep = "")

  SEASNames <- names(OutDF)[grep(pattern = "SEAS", names(OutDF))]
  FTNames <- names(OutDF)[grep(pattern = "FireType", names(OutDF))]
  SEAS_Matrix <- as.matrix(OutDF[, SEASNames])
  FT_matrix <- as.matrix(OutDF[, FTNames])
  Cols <- ncol(SEAS_Matrix)
  SEAS_Matrix[SEAS_Matrix == 0] <- NA
  Interval <- SEAS_Matrix[, 2:Cols] - SEAS_Matrix[, 1:Cols -
                                                    1]
  IntNames <- paste("INT", sprintf("%02d", 1:(Cols - 1)), sep = "")
  colnames(Interval) <- IntNames
  OutDF <- cbind(OutDF, Interval) %>% dplyr::select(-XYString)
  if (max_interval > 0) {
    FT_matrix <-
      fireTypeLowToHigh(
        max_interval = as.integer(max_interval),
        Interval_Matrix = Interval,
        Firetype_Matrix = FT_matrix
      )
    OutDF[, FTNames] <- FT_matrix
  }
  else if (max_interval < 0) {
    stop("max interfal cannot be less than 0")
  }
  else {

  }
  OutDF <- sf::st_as_sf(OutDF)
  LTR <- length(TimeSpan)
  SEAS_Matrix[is.na(SEAS_Matrix)] <- 0
  LBY <- matrix(NA, nrow(SEAS_Matrix), LTR)
  for (i in 1:LTR) {
    try({
      y = TimeSpan[i]
      LBY[, i] <- LBY_f(M = SEAS_Matrix, y)
    })
  }
  tYSF <- TimeSpan - t(LBY)
  YSF <- t(tYSF)
  YSFNames <- paste0("YSF", TimeSpan)
  LBYNames <- paste0("LBY", TimeSpan)
  LFTNames <- paste0("LFT", TimeSpan)
  colnames(YSF) <- YSFNames
  colnames(LBY) <- LBYNames
  print("calculating lookup matrix for getting last FireType by SEASON")
  SEAS_Matrix[SEAS_Matrix == 0] <- NA
  LUM <- matrix(NA, nrow(SEAS_Matrix), max.SEASON)
  for (i in 1:nrow(SEAS_Matrix)) {
    R <- i
    C <- as.numeric(stats::na.omit(SEAS_Matrix[i,]))
    V <- (FT_matrix[i, (1:length(C))])
    LUM[R, C] <- V
  }
  print("calculating last fire type")
  LFT <- matrix(NA, nrow(SEAS_Matrix), LTR)
  for (i in 1:nrow(SEAS_Matrix)) {
    LFT[i,] <- LUM[i, LBY[i,]]
  }
  colnames(LFT) <- LFTNames
  OutDF <- cbind(OutDF, YSF)
  OutDF <- cbind(OutDF, LBY)
  OutDF <- cbind(OutDF, LFT)
  OutDF$ID <- as.integer(rownames(OutDF))
  print("Completed making FH object")
  gc()
  results <-
    list(
      OutDF = OutDF,
      TimeSpan = TimeSpan,
      YSFNames = YSFNames,
      LBYNames = LBYNames,
      LFTNames = LFTNames
    )
  return(results)

}

#' FIRE History Processing for FAME stage1
#'
#' @details Takes a spatial polygon data file or r sf polygon dataset
#' checks that it contains the correct fields, and projection for FAME analysis
#' ( using helper function fhCheck) including fields FIRETYPE and SEASON
#'  Where polygons of different FIRETYPE or SEASON overlap the mainFHProcess()
#'  function  called in this function constructs unique
#'   non-overlapping polygon of their intersections ( and non intersecting areas
#'   ) and attributes each polygon with sequential fire SEASON (SEAS01, SEAS02
#'   ...) and corresponding FIRETYPE (TYPE01,TYPE02 ...)
#' @seealso [fhCheck()] for checking function run from within fhprocess1()
#'  [mainFHProcess]  main processing function to derive Fire history analysis for
#'  FAME
#'  [prepFH] helper function that splits larger inFH files into gridded subunits
#'   for parallel processing prior to passing to [mainFHProcess]
#'
#'
#' @param inFH Input fire history polygon/mutipolygon data set with columns
#'  SEASON and FIRETYPE  provided as shapefile, geopackage or ESRI geodatabase
#'  file, or the name of an sf  object. if the geometries contain other than
#'  polygon/mutipolygon an error will result.
#' @param inFHLayer Layer name if inFH has more than one layer
#' ( for instance in a .gpkg) this allows selection of a particular layer,
#' otherwise fist layer is used (Default = NULL)
#' @param OtherAndUnknown integer Value to use for cases where fire type is:
#'    "OTHER" or "UNKNOWN" = NA,
#'    "BURN" = 1,
#'    "BUSHFIRE" = 2.
#'    NA = Fire excluded from analysis.
#'  (usually set in settings file)
#' @param validFIRETYPE vector of valid names in the input FIRETYPE column in
#' the input fire history dataset(s), if the column contains NA or values not on
#' this list an error will occur
#' @param secondFH Second fire history to be combined with FH1 to make a fire
#'   scenario same formats as for inFH
#' @param secondFHLayer  Layer name if secondFH  has more than one layer ( for
#'   instance in a .gpkg) this allows selection of a particular layer, otherwise
#'   fist layer is used (Default = NULL)
#' @param baseFire Default NULL otherwise four digit integer SEASON  for fire
#'   applied #' across the whole bounding box
#' @param precision default 10000, precision value fed through to GEOS see
#'   [sf::st_set_precision], precision in metres = 1/precision 0 used precision
#'   of exiting files
#'
#' @return inFH1 an sf geometry collection of processed fire history for input
#'   into fhProcess2()
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @importFrom doSNOW registerDoSNOW
#' @export
#'
#' @examples
#' # randomFH<-generate_random_fire_history(20)
#' # plot(randomFH)
#' # outFH1<-fhProcess1(randomFH)
#' # plot(outFH1)
#' # #complete all fields in FH analysis
#' # outFH2<-fhProcess2(outFH1)
#' #
#' # plot(outFH2$OutDF,max.plot = 20)
#'
#'
fhProcess1<-function (inFH,
                      inFHLayer = NULL,
                      OtherAndUnknown = 2,
                      validFIRETYPE = c("BURN",
                                        "BUSHFIRE",
                                        "UNKNOWN",
                                        "OTHER"),
                      secondFH = NULL,
                      secondFHLayer = NULL,
                      baseFire = NULL,
                      precision =1)
{
  mySF <- fhCheck(inFH, inFHLayer, validFIRETYPE)
  if (!(is.null(secondFH) |length(secondFH)==0)) {
    mySF2 <- fhCheck(secondFH, secondFHLayer, validFIRETYPE)
    mySF <- dplyr::bind_rows(mySF, mySF2)
  }
  BBOX <- sf::st_cast(sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(mySF))),
                      "POLYGON")
  names(BBOX)[1] <- "geometry"
  if (is.null(baseFire)) {
  }
  else if (nchar(baseFire == 4)) {
    BBOX <- sf::st_as_sf(dplyr::bind_cols(BBOX, data.frame(SEASON = baseFire,
                                                           FIRETYPE = "BUSHFIRE")))

    mySF <- rbind(BBOX, mySF[, c("SEASON", "FIRETYPE")])

  }
  else {
    stop("baseFire argument is invalid")
  }
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BURN"] <- 1
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BUSHFIRE"] <- 2
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "OTHER"] <- OtherAndUnknown
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "UNKNOWN"] <- OtherAndUnknown

  st_precision(mySF)<- precision

  myPreppedFH <- prepFH(inFH =mySF)
  n_WithIntersects = sum(myPreppedFH$hasIntersects)
  n_WithoutIntersects = sum(!myPreppedFH$hasIntersects)
  print(paste("there are", n_WithIntersects, "gridcells with intersecting fires"))
  if (n_WithIntersects > 3) {
    cores<-parallel::detectCores()
    print(paste(cores, "cpu cores will be used in parallel"))
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)

    iterations<-n_WithIntersects
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    inFH1 <- foreach::foreach(i = myPreppedFH$hasData[myPreppedFH$hasIntersects],
                              .packages = c("sf", "dplyr", "tidyr", "magrittr",
                                            "FAMEFMR"),
                              .combine = dplyr::bind_rows,
                              .options.snow = opts) %dopar%
      {

        inFH1 <- mainFHProcess(x = i)
      }
  }
  else {
    print("parallel processing not being used")
    inFH1 <- foreach::foreach(i = myPreppedFH$hasData[myPreppedFH$hasIntersects],
                              .packages = c("sf", "dplyr", "tidyr", "magrittr",
                                            "FAMEFMR"), .combine = dplyr::bind_rows) %do%
      {
        inFH1 <- mainFHProcess(x = i)
      }
    return(inFH1)
  }
  if (n_WithoutIntersects > 0) {
    print("adding cells with no intersecting polygons")
    inFH1 <- dplyr::bind_rows(myPreppedFH$hasData[myPreppedFH$hasIntersects]) %>%
      dplyr::rename(SEAS01 = SEASON, FireType01 = FIRETYPE_NO) %>%
      dplyr::select(SEAS01, FireType01) %>% dplyr::bind_rows(inFH1)
  }
  return(inFH1)
}

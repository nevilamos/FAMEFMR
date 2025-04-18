#' FIRE History Processing for FAME stage1
#'
#' @details Takes a spatial polygon data file or r sf polygon collection
#' checks that it contains the correct fields, and projection for FAME analysis
#' ( using helper function checkFH) including fields FIRETYPE and SEASON
#'  Where polygons of different FIRETYPE or SEASON overlap the mainFHProcess()
#'  fucntion  called in this function constructs unique
#'   non-overlapping polygon of their intersections ( and non intersecting areas
#'   ) and attributes each polygon with sequential fire SEASON (SEAS01, SEAS02
#'   ...) and corresponding FIRETYPE (TYPE01,TYPE02 ...)
#' @seealso [checkFH()] for checking function run from within fhprocess1()
#'  [mainFHProcess]  main processing function to derive Fire history analysis for
#'  FAME
#'  [prepFH] helper function that splits larger inFH files into gridded subunits
#'   for parallel processing.
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
#' @param validFIRETYPE vector of valid firetype names in the input FIRETYPE
#' column, if the column contains NA or values not on this list an error will
#'  occur
#' @param secondFH Second fire history to be combined with FH1 to make a fire
#'  scenario same formats as for inFH
#' @param secondFHlayer  Layer name if secondFH  has more than one layer
#' ( for instance in a .gpkg) this allows selection of a particular layer,
#' otherwise fist layer is used (Default = NULL)
#' @param baseFire Default NULL otherwise four digit integer SEASON  for fire
#' applied #' across the whole bounding box
#'
#' @return inFH1 an sf geometry collection of processed fire history for input
#' into fhProcess2()
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
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
fhProcess1<-function(
  inFH,
  inFHLayer = NULL,
  OtherAndUnknown =2,
  validFIRETYPE = c("BURN", "BUSHFIRE", "UNKNOWN", "OTHER"),
  secondFH = NULL,
  secondFHLayer = NULL,
  baseFire = NULL){
  print(1)

  mySF<-fhCheck(inFH,inFHLayer,validFIRETYPE)
  if(!is.null(secondFH)){
    mySF2<-fhCheck(secondFH,secondFHLayer,validFIRETYPE)
    mySF<-dplyr::bind_rows(mySF,mySF2)
  }
  print(2)

  BBOX<-sf::st_cast(sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(mySF))),"POLYGON")
  names(BBOX)[1]<-"geometry"
  print(3)

  #Add fire everywhere in specified season
  if(is.null(baseFire)){
  }else if(nchar(baseFire==4)){
    #print(baseFire)
    BBOX<-sf::st_as_sf(dplyr::bind_cols(BBOX,data.frame(SEASON = baseFire,FIRETYPE="BUSHFIRE")))

    print(min(mySF$SEASON))
    mySF<-rbind(BBOX,mySF[,c("SEASON","FIRETYPE")])

    print(min(mySF$SEASON))
  } else{
    stop("baseFire argument is invalid")
  }
  print(4)
  # replacing FIRETYPE string names with integer values ----
  # could be optimised bu use of lookup with DT?
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BURN"] <- 1
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BUSHFIRE"] <- 2
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "OTHER"] <- OtherAndUnknown
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "UNKNOWN"] <- OtherAndUnknown
  print(5)



  myPreppedFH<-prepFH(inFH = mySF)
  #check the number of parts of the FH that have intersecting
  #polygons and chose whether to run in parallel and how many cores as a result
  n_WithIntersects = sum(myPreppedFH$hasIntersects)
  n_WithoutIntersects = sum(!myPreppedFH$hasIntersects)
  print(paste ("there are",n_WithIntersects,"gridcells with intersecting fires"))
  if(n_WithIntersects>3){
    cores <- parallel::detectCores()
    cores <- min(c(cores[1]-2,n_WithIntersects))
    print(paste(cores,"cpu cores will be used in parallel"))
    cl <- parallel::makeCluster(cores) #not to overload your computer
    doParallel::registerDoParallel(cl)
    try(on.exit(parallel::stopCluster(cl)))
    inFH1<-foreach::foreach(i=myPreppedFH$hasData[myPreppedFH$hasIntersects],.packages = c("sf","dplyr","tidyr","magrittr","FAMEFMR"),.combine=dplyr::bind_rows)%dopar%{
      print(i)

      inFH1<-mainFHProcess(x= i)
      # return(inFH1)
    }
  } else {
    print("parallel processing not being used")
    inFH1<-foreach::foreach(i=myPreppedFH$hasData[myPreppedFH$hasIntersects],.packages = c("sf","dplyr","tidyr","magrittr","FAMEFMR"),.combine=dplyr::bind_rows)%do%{

      inFH1<-mainFHProcess(x= i)
      #
    }
    return(inFH1)
  }

  #foreach loop to process each case where there are no intersections and bind into one  output SF
  if(n_WithoutIntersects>0){
    print("adding cells with no intersecting polygons")
    inFH1<-dplyr::bind_rows(myPreppedFH$hasData[myPreppedFH$hasIntersects]) %>%
      dplyr::rename(SEAS01=SEASON,FireType01=FIRETYPE_NO) %>%
      dplyr::select(SEAS01,FireType01) %>%
      dplyr::bind_rows(inFH1)
  }
  return(inFH1)
}

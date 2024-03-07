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
#' # plot(outFH2$outDF,max.plot = 20)
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
  mySF<-fhCheck(inFH,inFHLayer,validFIRETYPE)
  if(!is.null(secondFH)){
  mySF2<-fhCheck(secondFH,secondFHLayer)
  mySF<-dplyr::bind_rows(mySF,mySF2)
  }

  BBOX<-sf::st_cast(sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(mySF))),"POLYGON")
  names(BBOX)[1]<-"geometry"

  #Add fire everywhere in specified season
  if(is.null(baseFire)){
  }else if(nchar(baseFire==4)){
    BBOX<-dplyr::bind_cols(BBOX,data.frame(SEASON = baseFire,FIRETYPE="BUSHFIRE"))
    mySF<-dplyr::bind_rows(mySF,BBOX)
  } else{
    stop("baseFire argument is invalid")
  }
  # replacing FIRETYPE string names with integer values ----
  # could be optimised bu use of lookup with DT?
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BURN"] <- 1
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BUSHFIRE"] <- 2
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "OTHER"] <- OtherAndUnknown
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "UNKNOWN"] <- OtherAndUnknown



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




#' Main Fire History Fire Sequence analysis function
#'
#' @details The function takes sf polygon collection of fire history
#'   containing polygons with two fields: FIRETYPE and SEASON Where polygons of
#'   different FIRETYPE or SEASON overlap the function constructs unique
#'   non-overlapping polygon of their intersections ( and non intersecting areas
#'   ) and attributes each polygon with sequential fire SEASON (SEAS01, SEAS02
#'   ...) and corresponding FIRETYPE (TYPE01,TYPE02 ...)
#'
#' @param x input Fire History sf polygons with fields SEASON and FIRETYPE
#' @seealso [fhProcess1()] mainFHProcess() normally run from within it
#' @return output processed fire history with unique polygons each with sequence
#'  of fires and firetypes
#' @export
#'
#'
#' @examples
#' # randomFH<-generate_random_fire_history(20)
#' # plot(randomFH)
#' # #the function is run from within fhProcess1()
#' # outFH1<-fhProcess1(randomFH)
#' # plot(outFH1)
#' # #complete all fields in FH analysis
#' # outFH2<-fhProcess2(outFH1)
#' #
#' # plot(outFH2$outDF,max.plot = 20)
#'
#'
mainFHProcess<-function (x = inFH)  {
  mySF<-x %>%
    dplyr::select(SEASON,FIRETYPE_NO) %>%
    sf::st_make_valid(
      oriented = FALSE,
      #s2_options = s2::s2_options(snap = s2::s2_snap_precision(1e+07), ...),
      geos_method = "valid_structure",
      geos_keep_collapsed = FALSE
    )  %>%
    dplyr::arrange(SEASON,FIRETYPE_NO)

  SEASONS<-mySF$SEASON
  FIRETYPE_NOS<-mySF$FIRETYPE_NO

  outSF<-sf::st_intersection(mySF) %>%
    sf::st_collection_extract()%>%
    sf::st_cast(.,"POLYGON") %>%
    tidyr::unnest_wider(origins,names_sep = "_")




  Vnames<-names(outSF)[grep("origins",names(outSF))]

  originIndices<-outSF[,Vnames] %>%
    sf::st_drop_geometry() %>%
    as.matrix()

  seasonSequence<-as.data.frame(matrix(SEASONS[as.matrix(originIndices)],dim(originIndices)))

  paddedSequence<-formatC(1:dim(originIndices)[2],1,flag="0")

  names(seasonSequence)<-paste0("SEAS",paddedSequence)

  typeSequence<-as.data.frame(matrix(FIRETYPE_NOS[as.matrix(originIndices)],dim(originIndices)))

  names(typeSequence)<-paste0("FireType",paddedSequence)
  outSF<-dplyr::bind_cols(outSF,seasonSequence,typeSequence) %>%
    sf::st_as_sf() %>%
    dplyr::select(-all_of(c(Vnames,"SEASON","FIRETYPE_NO","n.overlaps")))%>%
    dplyr::group_by(across(c(-geometry))) %>%
    dplyr::summarise() %>%
    dplyr::ungroup()



  sf::st_crs(outSF)<-sf::st_crs(mySF)

  return(outSF)

}

#' Second stage of Fire history processing for FAME
#'
#' @description
#' Takes the results of the main FH process function  and calcutes
#'   Time Since fire (TSF) and Last Fire Type (LFT) and Last burnt year (LBY)
#'   for each SEASON as defined in the input arguments, these values are append
#'   to the output sf polygon dataframe.
#'
#'
#' @param inFH1 input partially processed Fire History from
#' @param start.SEASON first season for which output is wanted
#' (four digit year as integer), if NUll then second season in in history is used
#'  (cannot use first season because it has no interval, this may still fail if
#'   there is no overlap)
#' @param end.SEASON last season required, if NULL then maximum value in fire
#' history scenario used
#' @param max_interval maximum inter fire interval after a HIGH fire for which
#' subsequent fire is reassigned to HIGH, if 0 then no reassignment)
#'
#' @return A list containing: \itemize{
#' \item outDF A completed Fire Analysis Simple feature collection of polygons and/or
#'  multipolygons  with columns for Fire Season sequences, inter-fire intervals
#'  Last fire type,last fire SEASON and last burned SEASON for each season selected
#'  in the input range.
#' \item TimeSpan integer vector sequence of SEASONS to in the analysis output
#' \item YSFNames names of TSF years in output, needed by downstream functions
#' \item LBYNames names of LBY years in output, needed by downstream functions
#' \item LFTNames names of LBY years in output, needed by downstream functions }
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
#' # plot(outFH2$outDF,max.plot = 20)
#'
#'

fhProcess2<-function(inFH1,
                     start.SEASON = NULL,# first season for which output is wanted (four digit year as integer), if NULL then second season in in history is used (cannot use first season because it has no interval, this may still fail if there is no overlap)
                     end.SEASON = NULL,# last season required, if NULL then largest value in fire history scenario used
                     max_interval = 0 ){# maximum inter fire interval after a HIGH fire for which subsequent fire is reassigned to HIGH, if 0 then no reassignment)
  mySF<-inFH1

  mySF<-inFH1 %>%
    dplyr::group_by(across(-geometry)) %>%
    dplyr::summarise()
  # get just the names for the season and FireType sequence columns and make sure they are in ascending order
  SEASNames <- sort(names(mySF)[grep(pattern = "SEAS", names(mySF))])
  FTNames <- sort(names(mySF)[grep(pattern = "FireType", names(mySF))])

  # calculates all inter-fire intervals note by offsetting matrices by one column
  myDF<-mySF %>%
    sf::st_drop_geometry()

  SEAS_Matrix <- as.matrix(myDF[, SEASNames])

  FT_matrix <- as.matrix(myDF[, FTNames])

  Cols <- ncol(SEAS_Matrix)

  SEAS_Matrix[SEAS_Matrix == 0] <- NA

  Interval <- SEAS_Matrix[,2:Cols] - SEAS_Matrix[,1:Cols-1]

  IntNames <- paste("INT", sprintf("%02d", 1:(Cols-1)), sep = "")
  colnames(Interval) <- IntNames

  # bind the intervals calculated back to the main dataframe
  mySF <- cbind(mySF, Interval)[,c("geometry",SEASNames,FTNames,IntNames)]

  # work around for low firetype following after high firetype
  if(max_interval>0){
    FT_matrix<-fireTypeLowToHigh(max_interval = as.integer(max_interval),Interval_Matrix = Interval,Firetype_Matrix = FT_matrix)
    mySF[, FTNames]<-FT_matrix
  }else if (max_interval<0)  {
    stop ("max interval cannot be less than 0")
  }else {

  }





  mySF <- sf::st_as_sf(mySF)

  # time span (range of consecutive years) for which Fire History sequences are calculated
  # from start and end seasons defined in the settings file
  # else if NULL; start and end season for calculation is determined by the second smallest and max values in the dataset
  # the earliest start season is the date of a second season in the input data , next 3 lines prevent manual start
  # season setting less than this value, which would cause an error.
  uniqueSEASONS<-sort(unique(as.vector(SEAS_Matrix)))
  min.SEASON<-uniqueSEASONS[2] # second season in fire history
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
    max.SEASON <- max(uniqueSEASONS)
  } else {
    max.SEASON = end.SEASON
  }

  TimeSpan<-start.SEASON:max.SEASON



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

  # subtraction leads to matrix transposition which needs to be re transposed in the next couple of lines
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
  mySF <- cbind(mySF, YSF)
  mySF <- cbind(mySF, LBY)
  mySF <- cbind(mySF, LFT)
  mySF$ID <- as.integer(rownames(mySF))
  print("Completed making FH object")

  return(list("outDF" = mySF,"YSFNames" = YSFNames,"LBYNames" = LBYNames,"LFTNames" =LFTNames))
}




#' PrepareFire History
#'
#'@description
#' A function that handles the subdivision of a Fire History Dataset  for
#' Parallel processing in the FAME fire history processing workflow.
#' It calculates the number of grid cells to be used in parallelizing function.
#' splitting the input fire history into separate sf datasets on a grid with
#' areas of approximately 100 square km
#'
#' @param inFH
#'
#' @return A list containing: \itemize{
#' \item hasData list of one or more sf Fire History sub units of the input Fire
#' History containing fire boundaries, defined by a grid over the bounding box
#' of inFH with a grid of squreside*squareside cells
#' \item hasIntersects a logical vector indicating which elements of "hasData"
#' have intersecting Fire boundaries within them
#' \item squareSide integer number of units of one side of the grid used to divide
#' the inFH to produce "hasData" for parrallel compute.
#' Calculated as the floor() of area of inFH bounding box/10^8
#'  (ie a 10km square expressed in square metres)}
#' @export
#'
prepFH<-function(inFH = myFH){

  #calculate dimensions of grid sides
  squareside<-floor(sqrt(as.numeric(sf::st_area(sf::st_as_sfc(sf::st_bbox(inFH))))/10^8))
  squareside<-ifelse(squareside<1,1,squareside)
  if (squareside ==1){
    hasData = list(inFH)
  } else {
    grid<-sf::st_make_grid(inFH ,n=c(squareside,squareside))
    grid<-sf::st_as_sf(grid)
    grid_ID<-1:nrow(grid)
    grid$grid_ID<-grid_ID
    gridded<-sf::st_intersection(inFH,grid) %>% sf::st_make_valid()
    griddedList<-list()
    for(i in grid_ID){griddedList[[i]]<-gridded%>%
      dplyr::filter(grid_ID == i)}

    hasData<-griddedList[lapply(griddedList,nrow)>0]




  }
  hasIntersects<-unlist(lapply(hasData,FUN = function(x){
    max(unlist(lapply(as.list(sf::st_intersects(x)),length)))>0}))
  return(list("hasData"=hasData,
              "squareside" =squareside,
              "hasIntersects"=hasIntersects))
}

#' Basic Checks of input Fire History
#'
#' @description
#' Perfoms checks on input dataset to check it is correctly formatted for use in
#' FAME analysis
#'
#' @param inFH and input fire history polygon/mutipolygon data set with columns
#'  SEASON and FIRETYPE  provided as shapefile, geopackage or ESRI geodatabase
#'  file, or the name of an sf  object. if the geometries contain other than
#'  polygon/mutipolygon an error will result.
#' @param inFHLayer Layer name if inFH has more than one layer
#' ( for instance in a .gpkg) this allows selection of a particular layer,
#' otherwise first layer is used (Default = NULL)
#' #' @param validFIRETYPE vector of valid firetype names in the input FIRETYPE
#' column
#' @return correctly formatted Fire History Polygon dataset as sf for use in
#' fhProcess1()
#' @export
#'

fhCheck<-function(inFH=inFH,inFHLayer =inFHLayer,validFIRETYPE =validFIRETYPE){
  # Basic checks on input Fire History file----
  if("character" %in% (class(inFH))){
    if (tools::file_ext(inFH) %in% c("shp", "gdb", "gpkg")) {
      if (is.null(inFHLayer)) {
        mySF <- sf::st_read(dsn = inFH)
      }
      else {
        mySF <- sf::st_read(dsn = inFH, layer = inFHlayer)
      }
    }
  }   else if("sf" %in%  class(inFH)){
    mySF<-inFH} else {
      stop("inFH file is not a shapefile,geopackage or ESRI geodatabase nor is it a spatial features dataset")
    }

  if (!"sf" %in% class(mySF)) {
    stop("inFH file is not a spatial dataset")
  }
  if (!sf::st_geometry_type(mySF)[1] %in%
      c("POLYGON", "MULTIPOLYGON")) {
    stop("inFH file is not a polygon dataset")
  }
  mySFNames <- names(mySF)
  if ((!"SEASON" %in% mySFNames)){stop("inFH file does not contain field 'SEASON'")}

  if (!class(mySF[["SEASON"]]) %in% c("numeric", "integer")){stop("inFH file field 'SEASON' is not type integer")}
  if (anyNA(mySF$SEASON)) {stop("inFH file has missing values in  field 'SEASON'")}

  if ((!"FIRETYPE" %in% mySFNames))     {stop("inFH file does not contain field 'FIRETYPE'")}
  if (notAllIn(x = mySF$FIRETYPE, v = validFIRETYPE))  {stop("inFH shapefile has missing or invalid values in field 'FIRETYPE'")}

  inEPSG<-sf::st_crs(mySF)$epsg
  if (inEPSG == 3111){

  }else if(is.numeric(inEPSG)){
    mySF<-sf::st_transform(mySF,crs = "epsg:3111")
  } else{
    stop("input Fire History projection is not adequately defined. Ideally should be epsg:3111 (VicGrid94)")
  }
  return(mySF)
}




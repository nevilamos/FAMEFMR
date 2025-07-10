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
#' @param validFIRETYPE vector of valid names in the input FIRETYPE column in
#'   the input fire history dataset(s), if the column contains NA or values not
#'   on this list an error will occur
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
        mySF <- sf::st_read(dsn = inFH, layer = inFHLayer)
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
  sf::st_geometry(mySF)<-"geometry" # ensures that name of geometry column is always "geometry"
  return(mySF)
}


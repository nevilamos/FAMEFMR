#' also creates some rasters cropped to correct extent for instance for region and EFG
#' also gets indices of cells in raster of same extent as crop to the shape provided
#' @param REG_NO integer DELWP fire region number 1:6 ,99 for Statewide analysis,  or 7 for ad hoc boundary polygon default =7 (see look up table REG_LUT for values)
#' @param PUBLIC_LAND_ONLY Logical TRUE/FALSE
#' @param myRasterRes numeric raster resolution of the analysis in metres ( usually set in settings file or shiny app)
#' @param myPoly default clipPoly sf polygon data frame of LF_REGIONs (default) or ad hoc polygon - used in conjunction with REG_NO
#' @param generalRasterDir relative path to directory containing rasters of DELWP FIRE_REG, DELWP REGION, EFG, PUBLIC LAND (PLM_GEN)
#' @return A list containing:
#' \itemize{
#' \item Raster  raster cropped to extent of area of interest,
#' \item IDX integer vector cell numbers of cells in the cropped raster
#' \item clipIDX integer vector cell numbers only for cells with the input polygon
#' \item EFG integer vector EFG values for cells within clipped area
#' \item RGN integer vector Fire Region numbers for cells within clipped area
#' \item DELWP integer vector DELWP Region numbers for cells within clipped area
#' \item PLM logical for cells within clipped area
#' }
#' @export
cropNAborder  <- function(REG_NO = 7,
                          #see look up table REG_LUT for values
                          myRasterRes = RasterRes,
                          PUBLIC_LAND_ONLY,
                          myPoly = clipPoly,         #shapefile of LF_REGIONs (default) or adhoc polygon,
                          generalRasterDir = "./InputGeneralRasters"
){
  inputR <- inputRasters(myRasterRes)
  inR <- terra::rast(file.path(generalRasterDir, inputR[[1]]))
  Template <- inR
  terra::values(Template) <- 0
  idx<-cells(Template)


  #determines which file to use for masking to regions and values for cells in output region raster
  if(REG_NO %in% 1:6){
    Shape <- terra::vect(myPoly)
    Shape <- Shape[Shape$REGION_NO == REG_NO,]
    #cn <- tabularaster::cellnumbers(Template, Shape)
    inCells <- terra::cells(Template, Shape)[,"cell"]
    RGN <- Template
    terra::values(RGN)[cn] <- REG_NO
  } else  if(REG_NO == 99){
    Shape <- terra::vect(myPoly)
    inCells <- terra::cells(Template, Shape)[,"cell"]
    RGN <- inR
  } else  if(REG_NO == 7){
    Shape <- terra::vect(myPoly)
    inCells <- terra::cells(Template, Shape)[,"cell"]
    RGN <- Template
    terra::values(RGN)[inCells] <- REG_NO
  }

  # set parameters for extent of output
  Extent<-terra::ext(Shape)
  #  and crop rasters to output extent and  mask to selected area (RGN)
  RGNcrop <- terra::crop(RGN, Extent)
  #next line ensures extent aligned to RGN raster cells
  Extent<-terra::ext(RGNcrop)
  outCells<-terra::cells(RGNcrop,Shape)[,"cell"]

  inputRasters<-rast(file.path(generalRasterDir,unlist(inputR)))

  output<-as.list(as.data.frame(terra::values(inputRasters)[inCells,]))
  names(output)<-names(inputR)

  output$Raster<-rast(extent=Extent,res=225)
  output$inCells<-inCells
  output$outCells<-outCells



  # end of raster crop function
  return(output)
}



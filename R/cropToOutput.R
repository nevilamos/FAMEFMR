#'
#' Creates an output list containing vectors of input rasters, definintion ouf an empty
#' output raster (by extent and resoloution) cropped to the clipPolygon provided
#' also gets indices of cells in raster of same extent as crop to the shape provided
#' from input and output
#' @param REG_NO integer DELWP fire region number 1:6 ,99 for Statewide analysis,  or 7 for ad hoc boundary polygon default =7 (see look up table REG_LUT for values)
#' @param PUBLIC_LAND_ONLY Logical TRUE/FALSE
#' @param myRasterRes numeric raster resolution of the analysis in metres ( usually set in settings file or shiny app)
#' @param myPoly path of polygon file  LF_REGIONs (default) or ad hoc polygon - used in conjunction with REG_NO to crop the and output rasters.
#' @param generalRasterDir relative path to directory containing rasters of DELWP FIRE_REG, DELWP REGION, EFG, PUBLIC LAND (PLM_GEN)
#' @return A list containing:
#' \itemize{
#' \item Raster  raster cropped to extent of area of interest,
#' \item inCells integer vector cell numbers of cells in the inputRaster(s)
#'  that correspond to myPoly
#' \item outCells integer vector cell numbers of cells in the outputRaster(s)
#'  that correspond to the inputRaster
#' \item rasterDef an expression defining an empty raster with the correct extent
#' for the outputRaster
#' \item EFG integer vector EFG values for cells within clipped area
#' \item RGN integer vector Fire Region numbers for cells within clipped area
#' \item DELWP integer vector DELWP Region numbers for cells within clipped area
#' \item PLM logical for cells within clipped area
#' }
#' @export
cropToOutput  <- function(REG_NO = 7,
                          #see look up table REG_LUT for values
                          myRasterRes = RasterRes,
                          PUBLIC_LAND_ONLY,
                          myPoly = clipPoly,         #shapefile of LF_REGIONs (default) or adhoc polygon,
                          generalRasterDir = "./InputGeneralRasters"
){
  inputR <- inputRasters(myRasterRes)
  inR <- terra::rast(file.path(generalRasterDir, inputR[[1]]))
  Template <- inR
  terra::values(Template) <- NA


  #determines which file to use for masking to regions and values for cells in output region raster
  if(REG_NO %in% 1:6){
    Shape <- terra::vect(myPoly)
    Shape <- Shape[Shape$REGION_NO == REG_NO,]
    inCells <- terra::cells(Template, Shape)[,"cell"]
    RGN <- Template
    terra::values(RGN)[inCells] <- REG_NO
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

  # set  or extent of output
    #  and crop rasters to output extent and  mask to selected area (RGN)
  RGNTrim <- terra::trim(RGN)
  #next line ensures extent aligned to RGN raster cells
  Extent<-terra::ext(RGNTrim)

  outCells<-terra::cells(RGNTrim,Shape)[,"cell"]

  #get the values from all the standard input rasters cropped and masked to the region of interest masking to PLM if required
  inputRasters<-terra::rast(file.path(generalRasterDir,unlist(inputR)))
  #croppedInputs<-mask(crop(inputRasters,Extent),RGNTrim,maskvalues =c(NA,0))
  croppedInputs<-terra::crop(inputRasters,Extent)
  names(croppedInputs)<-names(inputR)
  #need to replace value for RGN which at the moment is take from the default regions
  IDX<-croppedInputs$IDX

  croppedInputs$RGN<-RGNTrim
  if (PUBLIC_LAND_ONLY ==TRUE){
    croppedInputs<-terra::mask(croppedInputs,croppedInputs$PLM)
    #need to reinstate all values of index not masked to public land only
    croppedInputs$IDX<-IDX
  }
  rm(IDX)



  output<-as.list(as.data.frame(terra::values(croppedInputs)))

  #output$Raster<-rast(extent=Extent,res=225)
  output$inCells<-inCells
  output$outCells<-outCells
  output$rasterDef<-parse(text =paste("terra::rast( extent = c(",Extent[1],",",Extent[2],",",Extent[3],",",Extent[4],")",",resolution =",myRasterRes,")"))
 # end of raster crop function
  return(output)
}



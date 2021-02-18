#' crop border of NA cells from rasters and get cell indices for remaning cells
#' function to get the minimum bounding box of the cells with non NA values in
#' a raster and save them to crop other rasters to same extent.
#' also creates some rasters cropped to correct extent for instance for region and EFG
#' also gets indices of cells in raster of same extent as crop to the shape provided
#' @param REG_NO integer DELWP fire region number 1:6 ,99 for Statewide analysis,  or 7 for ad hoc boundary polygon default =7 (see look up table REG_LUT for values)
#' @param PUBLIC_LAND_ONLY Logical TRUE/FALSE
#' @param myRasterRes numeric raster resolution of the analysis in metres ( usually set in settings file or shiny app)
#' @param myPoly default clipPoly sf polygon data frame of LF_REGIONs (default) or ad hoc polygon - used in conjunction with REG_NO
#' @param generalRasterDir relative path to directory containing rasters of DELWP FIRE_REG, DELWP REGION, EFG, PUBLIC LAND (PLM_GEN)
#' @return A list containing:
#' \itemize{
#' \item Raster  raster cropped of all border rows and columns that are all NA,
#' \item Extent extent of the raster
#' \item IDX integer vector cell numbers of cells in the cropped raster
#' \item clipIDX integer vector cell numbers only for cells with the input polygon
#' \item EFG integer vector EFG values for cells within clipped area
#' \item RGN integer vector Fire Region numbers for cells within clipped area
#' \item DELWP integer vector DELWP Region numbers for cells within clipped area
#' \item PLM logical for cells within clipped area
#' }
#' @export
cropNAborder  <- function(REG_NO = 7,              #see look up table REG_LUT for values
                          myRasterRes = RasterRes,
                          PUBLIC_LAND_ONLY,
                          myPoly = clipPoly,         #shapefile of LF_REGIONs (default) or adhoc polygon,
                          generalRasterDir = "./InputGeneralRasters"
){
  inputR <- inputRasters(myRasterRes)
  inR <- raster::raster(file.path(generalRasterDir, inputR$REGION.tif))
  Template <- inR
  raster::values(Template) <- NA

  #determines which file to use for masking to regions
  if(REG_NO %in% 1:6){
    Shape <- sf::read_sf(myPoly)
    Shape <- Shape[Shape$REGION_NO == REG_NO,]
    cn <- tabularaster::cellnumbers(Template, Shape)
    cn <- cn$cell_
    RGN <- Template
    raster::values(RGN)[cn] <- REG_NO
  }

  if(REG_NO == 99){
    Shape <- sf::read_sf(myPoly)
    cn <- tabularaster::cellnumbers(Template,Shape)
    cn <- cn$cell_
    RGN <- inR
    #raster::values(RGN)[cn] <- REG_NO
  }

  if(REG_NO == 7){
    Shape <- sf::read_sf(myPoly)
    cn <- tabularaster::cellnumbers(Template, Shape)
    cn <- cn$cell_
    RGN <- Template
    raster::values(RGN)[cn] <- REG_NO
  }
  # set parameters for
  x = RGN
  x.matrix <- is.na(raster::as.matrix(x))
  colNotNA <- which(colSums(x.matrix) != nrow(x))
  rowNotNA <- which(rowSums(x.matrix) != ncol(x))
  Extent <- raster::extent(x,
                   r1 = rowNotNA[1],
                   r2 = rowNotNA[length(rowNotNA)],
                   c1 = colNotNA[1],
                   c2 = colNotNA[length(colNotNA)]
  )

  # crop rasters
  RGN_ras <- raster::crop(RGN, Extent)
  FIREFMZ_ras <- raster::crop(raster::raster(file.path(generalRasterDir, inputR$FIREFMZ.tif)), Extent)
  DELWP_ras <- raster::crop( raster::raster(file.path(generalRasterDir, inputR$DELWP.tif)), Extent)
  EFG_ras <- raster::crop( raster::raster(file.path(generalRasterDir, inputR$EFG.tif)), Extent)
  IDX <- raster::values( raster::crop( raster::raster(file.path(generalRasterDir, inputR$IDX.tif)), Extent))
  PLM_ras <- raster::crop( raster::raster(file.path(generalRasterDir, inputR$PLM_GEN.tif)), Extent)
  EFG_ras <- raster::mask(EFG_ras, RGN_ras)

  # if choice has been made to restrict to public land (default) then EFG is masked to public land
  if(PUBLIC_LAND_ONLY == TRUE){
    EFG_ras <- raster::mask(EFG_ras, PLM_ras)
    DELWP_ras <- raster::mask(DELWP_ras, PLM_ras)
    FIREFMZ_ras <- raster::mask(FIREFMZ_ras, PLM_ras)
  }

  # set values as integers for rasters
  PLMVals <- as.integer(raster::values(PLM_ras))
  EFGvals <- as.integer(raster::values(EFG_ras))
  RGNvals <- as.integer(raster::values(RGN_ras))
  FIREFMZvals <- as.integer(raster::values(FIREFMZ_ras))
  DELWPvals <- as.integer(raster::values(DELWP_ras))

  # create list of values for function output
  output <- list("Raster" = RGN_ras,
                      "Extent" = Extent,
                      "clipIDX" = cn,
                      "EFG" = EFGvals,
                      "RGN" = RGNvals,
                      "IDX" = IDX,
                      "DELWP" = DELWPvals,
                      "FIREFMZ" = FIREFMZvals,
                      "PLM" = PLMVals
  )

  # end of raster crop function
  return(output)
}



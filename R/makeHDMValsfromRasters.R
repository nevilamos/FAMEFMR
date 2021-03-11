## FUNCTION makeHDMValsfromRasters  ------------------------------------------------------
#' Function makes a matrix of HDM values(1,NA) constrained to those cells that
#' are indexed in the cropped area
#' @param myHDMSpp_NO vector of TAXON_IDs for species to be included in output
#' @param myCropRasters list of rasters and indices and cell values created by function cropNAborder()
#' @importFrom iterators iter
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @export
makeHDMValsfromRasters <- function(myHDMSpp_NO = HDMSpp_NO,
                                   myCropRasters = cropRasters)
{
  HDMPaths <- dir(myCropRasters$HDM_RASTER_PATH, full.names=TRUE, pattern = ".tif$")
  HDMPaths <- HDMPaths[get_Spp_No(HDMPaths) %in%
                         myHDMSpp_NO]

  print("reading HDMvalues")

  # set up parallel processing
  cl<-parallel::makeCluster(Ncores)
  doParallel::registerDoParallel(cl, cores=Ncores)

  myHDMVals <- foreach::foreach(i = iterators::iter(HDMPaths),.combine = cbind,.packages = "raster") %dopar% {
    myVals <- raster::values(raster::raster(i))[myCropRasters$IDX]
    myVals}
  parallel::stopCluster(cl)

  colnames(myHDMVals) <- as.character(get_Spp_No(HDMPaths))
  return(myHDMVals)
}

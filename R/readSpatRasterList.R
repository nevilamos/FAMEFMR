
#' Read  qs file containing list with SpatRaster objects ( package terra)
#' @details Allows reloading of a list from and R session that contains SpatRasters the function recursively searches for PackedSpatRaster objects and coverts them to "SpatRaster" using terra::rast().  A work around to allow saving of sessions of R containing SpatRasters that can be reloaded ( saving the SpatRaster itself results in a null pointer error)
#' to save list containing SpatRasters from R use saveSpatRasterList
#' @param filePath file path to read the list object from ( as a package qs file with extension qs)
#'
#' @return list - if the input file contained saved PackedSpatRasters the list will contain the usable SpatRasters
#' @importFrom terra rast
#' @importFrom qs qread
#'
#' @export


readSpatRasterList<-function(filePath ="filepath to save list file to ending in qs"){
  stopifnot(tools::file_ext(filePath) == "qs")
  myList<-qs::qread(filePath)
  myList<-rapply(myList,
                 f=terra::rast,
                 classes = "PackedSpatRaster",
                 how="replace")
  return(myList)


}

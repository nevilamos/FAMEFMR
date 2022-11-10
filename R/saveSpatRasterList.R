
#' Save List containing SpatRaster objects ( package terra)
#' @details Allows saving and reloading of a list from and R session that contains SpatRasters the function recursively searches for SpatRaster objects and coverts them to "PackedSpatRaster" using terra::wrap()
#' A work around to allow saving of sessions of R containing SpatRasters that can be reloaded ( saving the SpatRaster itself results in a null pointer error)
# 'to read list back to R session use function readSpatRasterList
#' @param myList a list object that may contain SpatRasters
#' @param filePath file path to save the list object to ( as a package qs file with extention qs)
#'
#'
#' @importFrom terra wrap
#' @importFrom qs qsave
#'
#' @export
saveSpatRasterList<-function(myList,filePath ="filepath to save list file to ending in qs"){
  stopifnot(tools::file_ext(filePath) == "qs")

  myList<-rapply(myList,f=terra::wrap,classes = "SpatRaster",how="replace")
  qs::qsave(myList,filePath)

}


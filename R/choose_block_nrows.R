#' Choose a safe number of raster rows per block
#'
#' Computes a block height that respects both R's 2^31-1 vector limit
#' and a user-specified RAM budget.
#'
#' @param r SpatRaster
#' @param max_gb_block Maximum RAM (GB) per block
#' @param bytes_per_value Estimated bytes per value in memory
#'
#' @return Integer number of rows per block
#' @export
choose_block_nrows <- function(r, max_gb_block = 2, bytes_per_value = 4) {
  stopifnot(inherits(r, "SpatRaster"))

  nc <- ncol(r)
  nl <- nlyr(r)

  nrb_int <- floor(.Machine$integer.max / (nc * nl))
  nrb_int <- max(1L, as.integer(nrb_int))

  max_bytes <- max_gb_block * 1024^3
  max_elems <- floor(max_bytes / bytes_per_value)
  nrb_ram   <- floor(max_elems / (nc * nl))
  nrb_ram   <- max(1L, as.integer(nrb_ram))

  min(nrb_int, nrb_ram)
}

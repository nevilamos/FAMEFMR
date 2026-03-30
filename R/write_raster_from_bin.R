#' Write a packed binary representation back to a raster
#'
#' @param r_template SpatRaster template
#' @param bin_file Binary file from process_blocks_to_bin()
#' @param index_file Index CSV
#' @param max_ncol Number of layers in output
#' @param datatype Output GDAL datatype
#' @param NAflag NoData value
#' @param verbose Logical
#'
#' @export
write_raster_from_bin <- function(
    r_template,
    bin_file,
    index_file,
    #out_raster_file,
    max_ncol,
    datatype = "INT2U",
    NAflag   = 0,
    verbose  = TRUE
) {
  stopifnot(inherits(r_template, "SpatRaster"))
  stopifnot(file.exists(bin_file), file.exists(index_file))
  stopifnot(max_ncol >= 1)

  idx <- read.csv(index_file)
  nc <- ncol(r_template)

  #dir.create(dirname(out_raster_file), recursive = TRUE, showWarnings = FALSE)

  out_written_ok <- FALSE
  on.exit({
    if (!out_written_ok && file.exists(out_raster_file)) {
      unlink(out_raster_file, force = TRUE)
      unlink(paste0(out_raster_file, ".aux.xml"), force = TRUE)
    }
  }, add = TRUE)

  out <- terra::rast(r_template, nlyrs = as.integer(max_ncol))
  names(out) <- paste0("packed_", seq_len(max_ncol))
  f<-tempfile(fileext = ".tif")

  terra::writeStart(
    out,
    filename  = f,
    #overwrite = TRUE,
    wopt = list(
      datatype = as.character(datatype)[1],
      NAflag   = as.integer(NAflag)[1]
    )
  )

  con <- file(bin_file, "rb")
  on.exit(close(con), add = TRUE)

  for (k in seq_len(nrow(idx))) {
    rec <- idx[k, ]
    nrows_now <- rec$nrows
    cell_rows <- nrows_now * nc

    if (rec$ncol_out == 0L) {
      m_pad <- matrix(0L, cell_rows, max_ncol)
    } else {
      seek(con, rec$offset_bytes, origin = "start")
      v <- readBin(con, integer(), rec$n_values, size = 2,
                   endian = .Platform$endian, signed = TRUE)
      m <- matrix(v, cell_rows, rec$ncol_out, byrow = FALSE)
      m_pad <- matrix(0L, cell_rows, max_ncol)
      m_pad[, seq_len(rec$ncol_out)] <- m
    }

    terra::writeValues(
      out,
      v     = m_pad,
      start = rec$row_start,
      nrows = nrows_now
    )
  }

  terra::writeStop(out)
  out_written_ok <- TRUE

  out
}

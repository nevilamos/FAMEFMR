#' Process raster blocks into a packed binary representation
#'
#' Reads a SpatRaster block-wise, shifts non-zero values left per row,
#' drops trailing zero columns, and writes packed blocks to a binary file.
#'
#' @param r SpatRaster
#' @param out_prefix File prefix for output (.bin, .csv, .txt)
#' @param max_gb_block Maximum RAM per block (GB)
#' @param verbose Logical
#' @param progress Logical
#' @param progress_every Integer
#' @param show_rate_every Integer
#'
#' @details
#' Requires the Rcpp functions `shift_zero_in_place()` and
#' `drop_zero_cols()` to be available.
#'
#' @return Invisible list with paths and max_ncol
#' @export
process_blocks_to_bin <- function(
    r,
    out_prefix ="temp/packed_tmp",
    max_gb_block = 2,
    verbose = TRUE,
    progress = TRUE,
    progress_every = 1L,
    show_rate_every = 20L
) {
  stopifnot(inherits(r, "SpatRaster"))
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  nr <- nrow(r); nc <- ncol(r); nl <- nlyr(r)
  nrb <- choose_block_nrows(r, max_gb_block = max_gb_block, bytes_per_value = 12)

  bin_file   <- paste0(out_prefix, "_values_i16.bin")
  index_file <- paste0(out_prefix, "_index.csv")
  meta_file  <- paste0(out_prefix, "_meta.txt")

  con <- file(bin_file, open = "wb")
  on.exit(close(con), add = TRUE)

  idx <- data.frame(
    block = integer(),
    row_start = integer(),
    nrows = integer(),
    cell_rows = integer(),
    ncol_out = integer(),
    offset_bytes = numeric(),
    n_values = numeric(),
    stringsAsFactors = FALSE
  )

  max_ncol <- 0L
  n_blocks <- as.integer(ceiling(nr / nrb))

  pb <- NULL
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = n_blocks, style = 3)
    on.exit(try(utils::close(pb), silent = TRUE), add = TRUE)
  }

  terra::readStart(r)
  on.exit(terra::readStop(r), add = TRUE)

  row <- 1L
  block <- 1L

  while (row <= nr) {
    nrows_now <- min(nrb, nr - row + 1L)

    m <- terra::readValues(r, row = row, nrows = nrows_now, mat = TRUE)
    if (!is.integer(m)) m <- as.integer(m)
    if (anyNA(m)) m[is.na(m)] <- 0L

    shift_zero_in_place(m)
    m2 <- drop_zero_cols(m)

    ncol_out <- ncol(m2)
    max_ncol <- max(max_ncol, ncol_out)

    off <- seek(con)
    writeBin(as.integer(m2), con = con, size = 2, endian = .Platform$endian)

    cell_rows <- nrows_now * nc
    idx[nrow(idx) + 1L, ] <- list(
      block, row, nrows_now, cell_rows,
      ncol_out, off, cell_rows * ncol_out
    )

    if (progress && block %% progress_every == 0L)
      utils::setTxtProgressBar(pb, block)

    row <- row + nrows_now
    block <- block + 1L
  }

  write.csv(idx, index_file, row.names = FALSE)
  writeLines(
    c(
      paste0("nrow=", nr),
      paste0("ncol=", nc),
      paste0("nlyr_input=", nl),
      paste0("max_ncol_global=", max_ncol)
    ),
    meta_file
  )

  invisible(list(
    bin = bin_file,
    index = index_file,
    meta = meta_file,
    max_ncol = max_ncol
  ))
}

test_that("process_blocks_to_bin + write_raster_from_bin round-trip is correct", {

  skip_if_not_installed("terra")
  library(terra)

  # ---- mock required Rcpp functions ----
  # Shift: pack non-zeros left, preserve order (row-wise)
  shift_zero_in_place <- function(x) {
    nr <- nrow(x)
    nc <- ncol(x)
    for (i in seq_len(nr)) {
      row <- x[i, ]
      nz  <- row[row != 0L]
      if (length(nz) > 0) {
        x[i, ] <- c(nz, rep(0L, nc - length(nz)))
      } else {
        x[i, ] <- 0L
      }
    }
    invisible(NULL)
  }

  # Drop columns that are entirely zero
  drop_zero_cols <- function(x) {
    keep <- colSums(x != 0L) > 0L
    if (!any(keep)) {
      matrix(integer(0), nrow = nrow(x), ncol = 0)
    } else {
      x[, keep, drop = FALSE]
    }
  }

  # Inject mocks into process_blocks_to_bin environment
  env <- environment(process_blocks_to_bin)
  assign("shift_zero_in_place", shift_zero_in_place, envir = env)
  assign("drop_zero_cols", drop_zero_cols, envir = env)

  # ---- create a small test raster ----
  set.seed(42)
  r <- rast(nrows = 6, ncols = 5, nlyrs = 4)

  # Values: 0 = missing, positives are "events"
  vals <- sample(c(0L, 10L, 20L, 30L),
                 terra::ncell(r) * nlyrs(r),
                 replace = TRUE)

  values(r) <- vals

  td <- tempdir()
  out_prefix <- file.path(td, "packed_roundtrip")

  # ---- pass 1: process to binary ----
  res <- process_blocks_to_bin(
    r,
    out_prefix   = out_prefix,
    max_gb_block = 0.001,
    progress     = FALSE,
    verbose      = FALSE
  )

  expect_true(file.exists(res$bin))
  expect_true(file.exists(res$index))
  expect_true(res$max_ncol >= 0)

  # ---- pass 2: reconstruct raster ----
  out_file <- file.path(td, "packed_roundtrip.tif")

  write_raster_from_bin(
    r_template      = r,
    bin_file        = res$bin,
    index_file      = res$index,
    out_raster_file = out_file,
    max_ncol        = res$max_ncol,
    verbose = FALSE
  )

  expect_true(file.exists(out_file))

  out <- rast(out_file)

  # ---- structural checks ----
  expect_equal(nrow(out), nrow(r))
  expect_equal(ncol(out), ncol(r))
  expect_equal(nlyr(out), res$max_ncol)

  # ---- semantic round-trip check ----
  # For each cell:
  #   - take original values
  #   - drop zeros
  #   - compare with packed output (ignoring padded zeros)
  orig <- values(r, mat = TRUE)
  packed <- values(out, mat = TRUE)

  for (i in seq_len(nrow(orig))) {
    orig_row <- orig[i, ]
    nz <- orig_row[orig_row != 0L]

    if (length(nz) == 0L) {
      # packed row should be all zeros
      expect_true(all(packed[i, ] == 0L))
    } else {
      expect_equal(packed[i, seq_along(nz)], nz)
      expect_true(all(packed[i, (length(nz) + 1L):ncol(packed)] == 0L))
    }
  }

})

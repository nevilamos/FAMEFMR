test_that("process_blocks_to_bin creates consistent binary + index", {

  skip_if_not_installed("terra")

  library(terra)

  # ---- mock the required Rcpp functions ----
  # These mocks preserve the semantics needed for the test
  shift_zero_in_place <- function(x) {
    # no-op for testing
    invisible(NULL)
  }

  drop_zero_cols <- function(x) {
    # drop columns that are all zero
    keep <- colSums(x != 0L) > 0L
    if (!any(keep)) {
      matrix(integer(0), nrow = nrow(x), ncol = 0)
    } else {
      x[, keep, drop = FALSE]
    }
  }

  # Inject mocks into function environment
  env <- environment(process_blocks_to_bin)
  assign("shift_zero_in_place", shift_zero_in_place, envir = env)
  assign("drop_zero_cols", drop_zero_cols, envir = env)

  # ---- create a small test raster ----
  set.seed(1)
  r <- rast(nrows = 6, ncols = 5, nlyrs = 4)

  # values: zeros (missing) + small integers
  values(r) <- sample(c(0L, 10L, 20L),
                      terra::ncell(r) * nlyrs(r),
                      replace = TRUE)

  td <- tempdir()
  out_prefix <- file.path(td, "packed_test")

  # ---- run the function ----
  res <- process_blocks_to_bin(
    r,
    out_prefix   = out_prefix,
    max_gb_block = 0.001,
    progress     = FALSE,
    verbose      = FALSE
  )

  # ---- existence checks ----
  expect_true(file.exists(res$bin))
  expect_true(file.exists(res$index))
  expect_true(file.exists(res$meta))

  # ---- index structure checks ----
  idx <- read.csv(res$index)

  expect_true(all(c(
    "block", "row_start", "nrows",
    "cell_rows", "ncol_out",
    "offset_bytes", "n_values"
  ) %in% names(idx)))

  expect_true(all(idx$nrows > 0))
  expect_true(all(idx$cell_rows == idx$nrows * ncol(r)))
  expect_true(all(idx$ncol_out >= 0))
  expect_true(res$max_ncol >= max(idx$ncol_out))

  # ---- binary size consistency check ----
  # int16 => 2 bytes per value
  expected_bytes <- sum(idx$n_values) * 2
  actual_bytes <- file.info(res$bin)$size

  expect_equal(actual_bytes, expected_bytes)

})

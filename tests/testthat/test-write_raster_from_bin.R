test_that("write_raster_from_bin reconstructs raster", {
  r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 4)
  values(r) <- sample(c(0L, 10L, 20L), terra::ncell(r) * 4, replace = TRUE)

  td <- tempdir()
  out_prefix <- file.path(td, "packed_test")

  res <- process_blocks_to_bin(
    r,
    out_prefix = out_prefix,
    max_gb_block = 0.01,
    progress = FALSE,
    verbose = FALSE
  )

  out_file <- file.path(td, "packed_out.tif")

  write_raster_from_bin(
    r_template      = r,
    bin_file        = res$bin,
    index_file      = res$index,
    out_raster_file = out_file,
    max_ncol        = res$max_ncol,
    verbose = FALSE
  )

  out <- terra::rast(out_file)
  expect_equal(nlyr(out), res$max_ncol)
  expect_equal(nrow(out), n

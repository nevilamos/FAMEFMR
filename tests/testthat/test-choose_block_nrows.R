
test_that("choose_block_nrows returns positive integer", {
  r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 5)
  nrb <- choose_block_nrows(r, max_gb_block = 0.01)
  expect_true(is.integer(nrb))
  expect_gt(nrb, 0)
})

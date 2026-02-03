
test_that("make_epsg_numeric handles multiple formats", {

  expect_equal(make_epsg_numeric(3857), 3857)
  expect_equal(make_epsg_numeric("EPSG:3857"), 3857)
  expect_equal(make_epsg_numeric("3857"), 3857)
  expect_equal(make_epsg_numeric(sf::st_crs(3857)), 3857)

})



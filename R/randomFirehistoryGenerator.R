#' Generate Random Fire History Polygons
#'
#' @description Function generates random fire history polygon dataset
# for testing of
#'
#' @param nPolygons number of individual polygons to be generated
#' @param bbox bounding box within which the polygons are to be generated
#' including its crs is st_bbox format. default is a 20km square in Central
#' Eastern Victoria and epsg:3111 /VicGrid94 crs
#' @param maxPoints the maximum number of points used to define each polygon ( min value 11).
#' @param seasonRange vector of integer season values selected at random for polygon attribute SEASON
#' @param validFIRETYPE vector of valid names in the input FIRETYPE column in
#'   the input fire history dataset(s), if the column contains NA or values not
#'   on this list an error will occur
#' @param crs coordinate reference system  for the output set using sf_crs(),
#' default epsg:3111
#'
#' @return sf polygon
#' @export
#'
#' @examples
#' require(sf)
#' randomFH<-generate_random_fire_history(nPolygons=20,maxPoints=11)
#' plot(randomFH["SEASON"])
#' plot(randomFH["FIRETYPE"])
#'
generate_random_fire_history <- function(nPolygons =10,
                                         bbox = sf::st_bbox(c(
                                           xmin = 2600000,
                                           xmax = 2620000,
                                           ymax = 2420000,
                                           ymin = 2400000
                                         )),
                                         maxPoints = 11,
                                         seasonRange = 1980:2025,
                                         crs = sf::st_crs(3111),
                                         validFIRETYPE = c("BURN", "BUSHFIRE",
                                                           "UNKNOWN", "OTHER")
                                         ) {
  # Create an empty list to store the polygons
  polygon_list <- list()

  # Generate random polygons
  for (i in 1:nPolygons) {
    # Generate a random polygon with a random number of points
    xsize = bbox["xmax"] - bbox["xmin"]
    ysize = bbox["ymax"] - bbox["ymin"]
    nPoints <- sample(10:maxPoints, 1)

    polygon_coords <-
      cbind((runif(nPoints) * xsize * runif(1)) + runif(1) * xsize + bbox["xmin"],
            (runif(nPoints) * ysize * runif(1)) + runif(1) * ysize + bbox["ymin"]
      )



    # Add the polygon to the list
    polygon_list[[i]] <-
      sf::st_concave_hull(sf::st_multipoint(polygon_coords), .05)
  }
  sfc <- sf::st_sfc(polygon_list)

  # Convert sfc to sf object
  sf <- sf::st_sf(geometry = sfc)
  sf <- sf[sf::st_geometry_type(sf) == "POLYGON", ]
  sf[["SEASON"]] <- sample(seasonRange, nrow(sf), replace = T)
  sf[["FIRETYPE"]] <- sample(validFIRETYPE, nrow(sf), replace = T)
  sf <- sf %>% dplyr::arrange(SEASON, FIRETYPE)
  sf::st_crs(sf) <- crs

  return(sf)

}

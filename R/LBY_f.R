# Function to calculate last burnt year (LBY) from matrix of rows of fire season
# iterating by year (y) used in calc_TFI_2
#' Calculate last burned year matrix (LBY)
#' @details Function to calculate last burnt year (LBY) from matrix of rows of fire season
#'  iterating by year (y) used in calc_TFI_2
#' @param M numeric matrix of fire sequences sequence in rows, values are SEASON
#' @param y numeric  SEASON
#'
#' @return matrix  of last burned year row for each unique fire history column for each SEASON
#' @export
LBY_f <- function(M, y){
  M[M > y | M == 0] <- NA
  LBY <- Rfast::rowMaxs(M, value = TRUE)
  LBY[is.infinite(LBY)] <- NA
  return(LBY)
}

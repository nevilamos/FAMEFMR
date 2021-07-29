#' Calculate reformulated geomean objective
#'
#' @details Internal GSO function. Used within gso(), an input for nloptr()
#' calculates 'the reformulated objective'
#' @details The reformulated objective; nloptr expects a minimization problem.
#' The geomean objective, constraints (and corresponding gradients) all take the
#' same set of extra arguments.
#'
#' @param x a list/column, e.g. run$solution,
#' @param spp a data.frame, usually data or option data.frame
#'
#' @return returns a single value (the summed log of each row * x)
#' @export
#'

geomean.obj <- function (x, spp) {
  logsum <-
    sum(log(apply(spp, 1, function (row) {
      sum(x * row, na.rm = TRUE)
    })))
  return(-logsum)
}

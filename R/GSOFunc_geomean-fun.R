#' Calculate 'geomean'
#'
#' @details Internal GSO function. Used within Scenarios() and gso()
#' calculates 'geomean', the original objective. It is scaled by number of spp
#'
#' @details The geomean objective, constraints (and corresponding gradients) all take the
#' same set of extra arguments.
#'
#' @param x a list/column, e.g. run$solution,
#' @param spp a data.frame, usually data or option data.frame
#'
#' @return returns a single value (exponential of the summed log of each row * x,
#' divided by no. spp)
#' @export
#'

geomean.fun <- function(x, spp) {
  logsum <-
    sum(log(apply(spp, 1, function (row) {
      sum(x * row, na.rm = TRUE)
    })))
  return(exp(logsum / (dim(spp)[1])))
}

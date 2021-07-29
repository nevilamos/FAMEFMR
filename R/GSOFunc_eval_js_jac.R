#' Jacobian (partial derivative matrix) for the inequalities
#'
#' @details Internal GSO function. Used within gso(), an input for nloptr()
#'
#' @details Since each constraint is linear, the Jacobian is the coefficients of each
#' variable in eval_cs().  Calculates the coefficient for each variable in eval_cs()
#'
#'
#' @param x a list/column, e.g. run$solution,
#' @param spp a data.frame, usually data or option data.frame
#'
#' @return returns a vector
#' @export
#'

eval_cs_jac <- function(x, spp) {
  n <- length(StageNames)
  return(matrix(
    rep(c(rep(1, n), rep(0, length(
      x
    ))), length(x) / n),
    ncol = length(x),
    nrow = length(x) / n,
    byrow = TRUE
  ))
}

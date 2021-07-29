#' Calculate geomean objective gradient
#'
#' @details Internal GSO function. Used within gso(), an input for nloptr()
#' calculates the sum of the partial derivatives, with respect to each variable
#'
#' @details The reformulated objective; nloptr expects a minimization problem.
#' The geomean objective, constraints (and corresponding gradients) all take the
#' same set of extra arguments.
#'
#' @param x a list/column, e.g. run$solution,
#' @param spp a data.frame, usually data or option data.frame
#'
#' @return returns a single value (the sum of the partial derivatives)
#' @export
#'

geomean.grad <- function(x, spp) {
  ## The contribution of each row to the gradient.
  ## The derivative of log(sum(X*R)) wrt x is R[x]/(sum(X*R)).
  contrib <-
    apply(spp, 1, function(row) {
      -row / (sum(x * row, na.rm = TRUE))
    })
  ## Sum the contribution vectors to get the gradient
  grad <- apply(contrib, 1, function(dWt) {
    sum(dWt, na.rm = TRUE)
  })
  return(grad)
}

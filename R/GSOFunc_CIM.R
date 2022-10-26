#' Calculate mean and 95% credible interval
#'
#' @details Internal GSO function. Used within OptRunCI()
#' Handy function to calculate mean and 95% credible interval of input
#'
#' @param x a vector
#'
#' @return returns the mean of x, lower and upper bounds of the 95% credible interval of x
#' @export
#'

CIM <- function(x) {
  return(c(mean(x), stats::quantile(x, probs = c(0.025, 0.975))))
}

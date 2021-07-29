#' Inequality contraints replacement
#'
#' @details Internal GSO function. Used within gso(), an input for nloptr()
#' Replaces the proportion equality constraint with <=
#'
#' @details Since the objective is monotone, we can replace the proportion
#' equality constraint with <=.
#' Per-species limits should be do-able in the same form.
#'
#' @param x a list/column, e.g. run$solution,
#' @param spp a data.frame, usually data or option data.frame
#'
#' @return returns a vector
#' @export
#'

eval_cs <- function(x, spp) {
  n <- length(StageNames)
  Res <- rep(NA, length(x) / n)
  for (i in 1:length(Res)) {
    Res[i] <- sum(x[1:n + (i - 1) * n]) - 1
  }
  return(Res)
}

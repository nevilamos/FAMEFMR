#' Construct a n-length vector (c(k, ..., k))
#'
#' @details Internal GSO function. Used within props() and gso()
#' @details A handy helper function, constructs an n length vector, c(k, ... , k)
#'
#' @param n length of vector to create
#' @param k the value of the units within the n length vector
#'
#' @return returns a vector
#' @export
#'

const_vector <- function(n, k) {
  return(sapply(seq_len(n), function (x)
    return(k)))
}

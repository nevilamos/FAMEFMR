#' Shift zeros right in a matrix (0 treated as empty)
#'
#' Treats \code{0L} as "empty" and compacts non-zero entries along the algorithm's
#' traversal order, filling remaining positions with \code{0L}.
#'
#' @param x A matrix coercible to integer. Zeros (\code{0L}) are treated as empty.
#' @return An integer matrix with zeros shifted.
#' @export
shift_zero <- function(x) {
  if (!is.matrix(x)) {
    stop("`x` must be a matrix.", call. = FALSE)
  }
  if (!is.integer(x)) {
    x <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x),
                dimnames = dimnames(x))
  }
  shift_zero_cpp(x)
}


#' Shift zeros in place (C++ implementation)
#'
#' Treats \code{0L} as "empty" and compacts non-zero entries along the algorithm's
#' traversal order, filling remaining positions with \code{0L}.
#'
#' @details
#' This function calls a C++ implementation via Rcpp. It is intended for
#' performance-critical workflows.
#'
#' **In-place semantics:** R uses copy-on-modify. This function *attempts* to
#' modify \code{x} in place, but R may duplicate \code{x} if it is referenced
#' elsewhere. For guaranteed behavior, use \code{\link{shift_zero}} which returns
#' a new matrix.
#'
#' @param x A matrix coercible to integer. Zeros (\code{0L}) are treated as empty.
#'
#' @return Invisibly returns \code{x} (an integer matrix) after shifting zeros.
#' The primary effect is the side-effect on \code{x}.
#'
#' @seealso \code{\link{shift_zero}} for a pure functional (copying) version.
#'
#' @export
shift_zero_in_place <- function(x) {
  if (!is.matrix(x)) {
    stop("`x` must be a matrix.", call. = FALSE)
  }

  # Coerce to integer matrix if needed, preserving dimnames
  if (!is.integer(x)) {
    x <- matrix(as.integer(x),
                nrow = nrow(x), ncol = ncol(x),
                dimnames = dimnames(x))
  }

  shift_zero_in_place_cpp(x)

  invisible(x)
}


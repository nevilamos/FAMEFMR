#' Drop all-zero columns (C++ implementation)
#'
#' Removes columns that contain only zeros (\code{0L}) from an integer matrix.
#'
#' @details
#' This function calls a C++ implementation via Rcpp. It is intended for
#' performance-critical workflows.
#'
#' By default, a column is retained if it contains **any** non-zero entry.
#'
#' @param mat A matrix coercible to integer.
#'
#' @return An integer matrix containing only columns with at least one non-zero
#' value.
#'
#' @export
drop_zero_cols <- function(mat) {
  if (!is.matrix(mat)) {
    stop("`mat` must be a matrix.", call. = FALSE)
  }

  # Coerce to integer matrix if needed, preserving dimnames
  if (!is.integer(mat)) {
    mat <- matrix(as.integer(mat),
                  nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = dimnames(mat))
  }

  drop_zero_cols_cpp(mat)
}


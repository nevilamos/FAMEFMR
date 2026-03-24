#' Split fixed-width integers into leading and trailing digit groups
#'
#' Splits each element of an integer vector or integer matrix into two parts:
#' the first \code{k_first} leading digits and the last \code{k_last} trailing digits,
#' assuming each value has exactly \code{n_digits} decimal digits (based on absolute value).
#'
#' A strict digit-length check is applied to every non-NA element. If any element
#' does not have exactly \code{n_digits} digits, the function errors and reports
#' the first offending position (index for vectors; row/col for matrices).
#'
#' Trailing groups that include leading zeros are returned as integers, so those
#' leading zeros are naturally dropped. For example, the last 2 digits of
#' \code{123001} are returned as \code{1} (i.e., "01" becomes 1).
#'
#' @param x An integer vector or integer matrix of fixed-width identifiers.
#'   Must be of type \code{integer}; numeric/double inputs are rejected to avoid truncation.
#' @param n_digits Integer scalar. Total number of digits each element of \code{x} must have.
#' @param k_first Integer scalar. Number of leading digits to extract.
#' @param k_last Integer scalar. Number of trailing digits to extract.
#'
#' @return A list with components \code{first} and \code{last}.
#' If \code{x} is a vector, both are integer vectors of the same length.
#' If \code{x} is a matrix, both are integer matrices with the same dimensions as \code{x}.
#'
#' @examples
#' # Vector input
#' x <- c(123001L, 999900L)
#' split_integer_first_last(x, n_digits = 6, k_first = 4, k_last = 2)
#'
#' # Matrix input
#' m <- matrix(c(123001L, 999900L,
#'               123099L, 999901L), nrow = 2, byrow = TRUE)
#' split_integer_first_last(m, n_digits = 6, k_first = 4, k_last = 2)
#'
#' # Numeric input rejected (avoids truncation):
#' \dontrun{
#' split_integer_first_last(c(123001, 999900), n_digits = 6, k_first = 4, k_last = 2)
#' }
#'
#' @export
split_integer_first_last <- function(x, n_digits, k_first, k_last) {
  if (!is.integer(x)) {
    stop(
      "`x` must be an integer vector or integer matrix (use the `L` suffix, e.g. 123001L). ",
      "Numeric/double inputs are not accepted to avoid truncation.",
      call. = FALSE
    )
  }

  n_digits <- as.integer(n_digits)
  k_first  <- as.integer(k_first)
  k_last   <- as.integer(k_last)

  if (length(n_digits) != 1L || is.na(n_digits)) stop("`n_digits` must be a single non-NA integer.", call. = FALSE)
  if (length(k_first)  != 1L || is.na(k_first))  stop("`k_first` must be a single non-NA integer.", call. = FALSE)
  if (length(k_last)   != 1L || is.na(k_last))   stop("`k_last` must be a single non-NA integer.", call. = FALSE)

  if (is.matrix(x)) {
    split_integer_first_last_matrix_cpp(x, n_digits = n_digits, k_first = k_first, k_last = k_last)
  } else {
    split_integer_first_last_cpp(x, n_digits = n_digits, k_first = k_first, k_last = k_last)
  }
}

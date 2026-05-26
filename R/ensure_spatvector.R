#' Ensure input is a SpatVector
#'
#' Accepts either a `terra::SpatVector` object or a file path that can be
#' read into a `SpatVector`. Returns a valid `SpatVector`, or stops with a
#' warning message if unsuccessful.
#'
#' @param x A `terra::SpatVector` object or a character file path.
#'
#' @return A `terra::SpatVector`.
#'
#' @examples
#' library(terra)
#'
#' v <- vect(system.file("ex/lux.shp", package = "terra"))
#' ensure_spatvector(v)
#'
#' ensure_spatvector(system.file("ex/lux.shp", package = "terra"))
#'
#' @export
ensure_spatvector <- function(x) {

  # Case 1: already a SpatVector
  if (inherits(x, "SpatVector")) {
    return(x)
  }

  # Case 2: character path
  if (is.character(x) && length(x) == 1) {

    if (!file.exists(x)) {
      warning(sprintf("File does not exist: %s", x))
      stop("Input could not be converted to SpatVector.", call. = FALSE)
    }

    v <- tryCatch(
      terra::vect(x),
      error = function(e) {
        warning(sprintf("Failed to read file as SpatVector: %s", e$message))
        return(NULL)
      }
    )

    if (inherits(v, "SpatVector")) {
      return(v)
    }
  }

  # Fallback: invalid input
  warning("Input is neither a SpatVector nor a readable file path.")
  stop("Input could not be converted to SpatVector.", call. = FALSE)
}

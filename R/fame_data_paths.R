#' Return the root FAME data directory
#'
#' The directory is obtained first from the `FAME.data_dir`
#' option and then from the `FAME_DATA_DIR` environment
#' variable.
#'
#' @param must_work Logical. Should the directory already exist?
#'
#' @return Normalised path to the external FAME data directory.
#' @export
fame_data_dir <- function(must_work = TRUE) {

  data_dir <- getOption(
    "FAME.data_dir",
    Sys.getenv("FAME_DATA_DIR", unset = "")
  )

  if (!is.character(data_dir) ||
      length(data_dir) != 1L ||
      is.na(data_dir) ||
      !nzchar(data_dir)) {

    stop(
      paste0(
        "The FAME external data directory has not been configured.\n",
        "Set it with either:\n\n",
        "  Sys.setenv(FAME_DATA_DIR = 'C:/path/to/FAME_data')\n\n",
        "or:\n\n",
        "  options(FAME.data_dir = 'C:/path/to/FAME_data')"
      ),
      call. = FALSE
    )
  }

  data_dir <- path.expand(data_dir)

  if (must_work && !dir.exists(data_dir)) {
    stop(
      "FAME_DATA_DIR does not exist: ",
      data_dir,
      call. = FALSE
    )
  }

  normalizePath(
    data_dir,
    winslash = "/",
    mustWork = must_work
  )
}


#' Construct a path within the FAME data directory
#'
#' @param ... Directory names or filenames below FAME_DATA_DIR.
#' @param must_work Logical. Should the root data directory exist?
#'
#' @return A platform-independent file path.
#' @export
fame_data_path <- function(..., must_work = TRUE) {

  file.path(
    fame_data_dir(must_work = must_work),
    ...
  )
}

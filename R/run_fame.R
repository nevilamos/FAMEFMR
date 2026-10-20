#' Run the FAME Shiny application
#'
#' Launches the FAME Shiny application bundled with FAMEFMR.
#'
#' @param data_dir Optional directory containing writable FAME data.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return The value returned invisibly by [shiny::runApp()].
#'
#' @export
run_fame <- function(
    data_dir = getOption(
      "FAME.data_dir",
      Sys.getenv("FAME_DATA_DIR", unset = "")
    ),
    ...
) {

  app_dir <- system.file(
    "shinyapp",
    package = "FAMEFMR"
  )

  if (!nzchar(app_dir)) {
    stop(
      "FAME Shiny application could not be found.",
      call. = FALSE
    )
  }

  if (nzchar(data_dir)) {

    data_dir <- normalizePath(
      data_dir,
      mustWork = FALSE
    )

    dir.create(
      data_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    options(FAME.data_dir = data_dir)

    Sys.setenv(
      FAME_DATA_DIR = data_dir
    )
  }

  shiny::runApp(
    appDir = app_dir,
    ...
  )
}

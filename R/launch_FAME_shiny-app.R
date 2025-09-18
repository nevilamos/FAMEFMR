#' Launch FAME shiny app
#' @export
launchFAMEShiny <- function() {
  appDir <- system.file("application", package = "FAMEFMR")
  if (appDir == "") {
    stop("Could not find FAMEShiny. Try re-installing `FAMEFMR`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

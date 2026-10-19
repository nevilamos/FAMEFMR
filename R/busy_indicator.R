#' Busy indicator UI wrapper for Shiny action buttons
#'
#' Wrap a Shiny button with a loading indicator, success checkmark,
#' and error display for improved user experience.
#'
#' @param button A Shiny button tag (e.g. from \code{actionButton()}).
#'
#' @return A \code{shiny.tag} object containing the modified button UI.
#'
#' @details
#' This function must be used together with \code{withBusyIndicatorServer()}.
#' It injects loading, success, and error indicators linked to the button ID.
#'
#' @examples
#' \dontrun{
#' withBusyIndicatorUI(actionButton("go", "Run"))
#' }
#'
#' @import shiny
#' @importFrom shinyjs hidden
#' @export
withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      hidden(
        img(src = "ajax-loader.gif", class = "btn-loading-indicator"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("circle-exclamation"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

#' Busy indicator server wrapper for Shiny actions
#'
#' Execute an expression while displaying loading, success, and error states
#' for a corresponding button created with \code{withBusyIndicatorUI()}.
#'
#' @param buttonId Character. The ID of the button.
#' @param expr Expression to evaluate when the button is triggered.
#'
#' @return The result of evaluating \code{expr}, or triggers error handling.
#'
#' @details
#' Handles UI updates via \pkg{shinyjs}, including:
#' \itemize{
#'   \item Showing a loading spinner
#'   \item Displaying a success checkmark briefly
#'   \item Showing an error message if evaluation fails
#' }
#'
#' @examples
#' \dontrun{
#' observeEvent(input$go, {
#'   withBusyIndicatorServer("go", {
#'     Sys.sleep(2)
#'     "done"
#'   })
#' })
#' }
#'
#' @import shiny
#' @importFrom shinyjs show hide enable disable delay
#' @export
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  disable(buttonId)
  show(selector = loadingEl)
  hide(selector = doneEl)
  hide(selector = errEl)
  on.exit({
    enable(buttonId)
    hide(selector = loadingEl)
  })

  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    show(selector = doneEl)
    delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                              time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

#' Handle errors for busy indicator buttons
#'
#' Display an error message associated with a button when an error occurs
#' during execution of \code{withBusyIndicatorServer()}.
#'
#' @param err An error object.
#' @param buttonId Character. The ID of the button.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @importFrom shinyjs html show
#' @export
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}

#' Default CSS for busy indicator buttons
#'
#' CSS styles required for proper display of the busy indicator UI elements.
#'
#' @return A character string containing CSS.
#'
#' @examples
#' \dontrun{
#' tags$style(appCSS)
#' }
#'
#' @export
appCSS <- "
.btn-loading-container {
  margin-left: 10px;
  font-size: 1.2em;
}
.btn-done-indicator {
  color: green;
}
.btn-err {
  margin-top: 10px;
  color: red;
}
"

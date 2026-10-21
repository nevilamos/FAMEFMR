#' Add a busy indicator to a Shiny button
#'
#' Wraps a Shiny button with UI elements used to display a loading indicator,
#' success indicator, and error message while a server-side operation is
#' running.
#'
#' This function should be used together with
#' [withBusyIndicatorServer()] on the server side.
#'
#' The wrapped button must have an `id`, as this is used to associate the
#' loading, completion, and error indicators with the corresponding button.
#'
#' @param button A Shiny button UI element, typically created using
#'   [shiny::actionButton()] or a similar function.
#'
#' @return A Shiny `div` containing the supplied button together with the
#'   loading, completion, and error indicator elements.
#'
#' @details
#' The function expects an image named `ajax-loader.gif` to be available to
#' the Shiny application as a web-accessible resource.
#'
#' The associated CSS defined in `appCSS` should also be included in the
#' application's UI.
#'
#' This implementation is based on the busy indicator example from Dean Attali's
#' Advanced Shiny examples:
#' \url{https://github.com/daattali/advanced-shiny/blob/master/busy-indicator/helpers.R}
#'
#' @seealso [withBusyIndicatorServer()], [errorFunc()]
#'
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


#' Run an expression with a button busy indicator
#'
#' Runs an expression while displaying a loading indicator associated with a
#' Shiny button. The button is temporarily disabled while the expression runs.
#'
#' If the expression completes successfully, a completion indicator is shown
#' briefly. If an error occurs, the error is passed to [errorFunc()] and
#' displayed alongside the button.
#'
#' This function is the server-side counterpart to
#' [withBusyIndicatorUI()].
#'
#' @param buttonId Character string giving the input ID of the Shiny button.
#'
#' @param expr An expression to evaluate while the busy indicator is active.
#'
#' @return The value returned by `expr` if it completes successfully.
#'   If an error occurs, the error is handled by [errorFunc()].
#'
#' @details
#' When called, the function:
#'
#' \enumerate{
#'   \item Disables the button.
#'   \item Displays the loading indicator.
#'   \item Evaluates `expr`.
#'   \item Displays a completion indicator if evaluation succeeds.
#'   \item Displays an error message if evaluation fails.
#'   \item Re-enables the button when processing finishes.
#' }
#'
#' The function relies on functions from the `shinyjs` package for controlling
#' UI elements.
#'
#' This implementation is based on the busy indicator example from Dean Attali's
#' Advanced Shiny examples:
#' \url{https://github.com/daattali/advanced-shiny/blob/master/busy-indicator/helpers.R}
#'
#' @seealso [withBusyIndicatorUI()], [errorFunc()]
#'
#' @export
withBusyIndicatorServer <- function(buttonId, expr) {

  # UX stuff: show the "busy" message, hide the other messages,
  # disable the button
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

  # Try to run the code when the button is clicked and show an error
  # message if an error occurs or a success message if it completes
  tryCatch({
    value <- expr

    show(selector = doneEl)

    delay(
      2000,
      shinyjs::hide(
        selector = doneEl,
        anim = TRUE,
        animType = "fade",
        time = 0.5
      )
    )

    value

  }, error = function(err) {
    errorFunc(err, buttonId)
  })
}


#' Display an error associated with a busy-indicator button
#'
#' Displays an error message in the error container associated with a button
#' created using [withBusyIndicatorUI()].
#'
#' This function is primarily an internal helper used by
#' [withBusyIndicatorServer()] when an error occurs while evaluating the
#' supplied expression.
#'
#' @param err An error condition object. The error message is obtained from
#'   `err$message`.
#'
#' @param buttonId Character string giving the input ID of the Shiny button
#'   associated with the error.
#'
#' @return Called for its side effects. The function updates and displays
#'   the error element associated with the button.
#'
#' @details
#' The prefix `"ddpcr: "` is removed from the error message, if present,
#' before the message is displayed.
#'
#' The function uses [shinyjs::html()] and [shinyjs::show()] to update and
#' display the error message.
#'
#' This implementation is based on the busy indicator example from Dean Attali's
#' Advanced Shiny examples:
#' \url{https://github.com/daattali/advanced-shiny/blob/master/busy-indicator/helpers.R}
#'
#' @seealso [withBusyIndicatorUI()], [withBusyIndicatorServer()]
#'
#' @keywords internal
errorFunc <- function(err, buttonId) {

  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)

  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)

  shinyjs::html(
    html = errMessage,
    selector = errElMsg
  )

  shinyjs::show(
    selector = errEl,
    anim = TRUE,
    animType = "fade"
  )
}


#' CSS for Shiny busy indicators
#'
#' CSS definitions used by [withBusyIndicatorUI()] to style loading,
#' completion, and error indicators.
#'
#' The CSS should be included in the Shiny application's UI, for example using
#' `shiny::tags$style(appCSS)`.
#'
#' @format A character string containing CSS definitions.
#'
#' @seealso [withBusyIndicatorUI()], [withBusyIndicatorServer()]
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

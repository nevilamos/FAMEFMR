#' Add a busy indicator to a Shiny button
#'
#' Wraps a Shiny button with UI elements used to display a loading indicator,
#' success indicator, and error message while a server-side operation is
#' running.
#'
#' This function should be used together with [withBusyIndicatorServer()] on
#' the server side.
#'
#' The supplied button must have an \code{id} attribute. This ID is used to
#' associate the loading, completion, and error indicators with the
#' corresponding button.
#'
#' @param button A Shiny button UI element, typically created with
#'   [shiny::actionButton()] or a similar function. The button must have an
#'   \code{id} attribute.
#'
#' @return A Shiny \code{div} containing the supplied button together with
#'   loading, completion, and error indicator elements.
#'
#' @details
#' The function expects an image named \code{ajax-loader.gif} to be available
#' to the Shiny application as a web-accessible resource.
#'
#' The CSS definitions in [appCSS] should also be included in the application's
#' UI.
#'
#' This implementation is based on the busy-indicator example from Dean
#' Attali's Advanced Shiny examples:
#' \url{https://github.com/daattali/advanced-shiny/blob/master/busy-indicator/helpers.R}
#'
#' @seealso [withBusyIndicatorServer()], [errorFunc()], [appCSS]
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
      div(
        class = "btn-err",
        div(
          icon("circle-exclamation"),
          tags$b("Error: "),
          span(class = "btn-err-msg")
        )
      )
    )
  )
}


#' Run an expression with a button busy indicator
#'
#' Evaluates an expression while displaying a loading indicator associated
#' with a Shiny button. The button is temporarily disabled while the
#' expression is evaluated.
#'
#' If the expression completes successfully, a completion indicator is shown
#' briefly. If an error occurs, the error is passed to [errorFunc()] and
#' displayed alongside the button.
#'
#' This function is the server-side counterpart to [withBusyIndicatorUI()].
#'
#' @param buttonId A character string giving the input ID of the Shiny button.
#'
#' @param expr An expression to evaluate while the busy indicator is active.
#'
#' @return The value returned by \code{expr} if evaluation completes
#'   successfully. If an error occurs, the error is handled by [errorFunc()]
#'   and the value returned by the error handler is returned.
#'
#' @details
#' When called, the function:
#'
#' \enumerate{
#'   \item Disables the button.
#'   \item Displays the loading indicator.
#'   \item Evaluates \code{expr}.
#'   \item Displays a completion indicator if evaluation succeeds.
#'   \item Displays an error message if evaluation fails.
#'   \item Re-enables the button when processing finishes.
#' }
#'
#' The function relies on functions from the \pkg{shinyjs} package to control
#' the button and associated indicator elements. [withBusyIndicatorUI()] must
#' be used for the corresponding button so that the required UI elements and
#' \code{data-for-btn} attribute are present.
#'
#' This implementation is based on the busy-indicator example from Dean
#' Attali's Advanced Shiny examples:
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
#' This function is an internal helper used by [withBusyIndicatorServer()]
#' when an error occurs while evaluating the supplied expression.
#'
#' @param err An error condition object. The error message is obtained from
#'   \code{err$message}.
#'
#' @param buttonId A character string giving the input ID of the Shiny button
#'   associated with the error.
#'
#' @return Called for its side effects. Returns the result of
#'   [shinyjs::show()] invisibly while updating and displaying the error
#'   element associated with the button.
#'
#' @details
#' The prefix \code{"ddpcr: "}, if present at the beginning of the error
#' message, is removed before the message is displayed.
#'
#' The function uses [shinyjs::html()] to insert the error message and
#' [shinyjs::show()] to display the associated error container.
#'
#' This implementation is based on the busy-indicator example from Dean
#' Attali's Advanced Shiny examples:
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
#' CSS definitions used by [withBusyIndicatorUI()] to style the loading,
#' completion, and error indicators.
#'
#' The CSS should be included in the Shiny application's UI, for example with
#' \code{shiny::tags$style(appCSS)}.
#'
#' @format A character string containing CSS definitions for the busy
#'   indicator UI elements.
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

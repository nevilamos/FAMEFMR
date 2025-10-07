#' Validated Numeric Input Module - UI
#'
#' Provides a numeric input that validates user input on Enter or blur.
#' Accepts NA (if allowed) or numeric/integer values within optional bounds.
#' Invalid input triggers a small popup modal and resets to the last valid value.
#'
#' Label and validation messages automatically adapt based on `integerOnly`, `minVal`, and `maxVal`.
#'
#' @param id Module ID.
#' @param integerOnly Logical. If TRUE (default), label and validation refer to "integer"; otherwise "number".
#' @param label Optional label text to display. If NULL, label is automatically generated using minVal and maxVal.
#' @param minVal Optional minimum allowed value (used for display and validation).
#' @param maxVal Optional maximum allowed value (used for display and validation).
#' @param defaultVal Optional default value (used as starting value in the input, defaults to NA).
#' @return Shiny UI element for the numeric input.
#' @examples
#' \dontrun{
#' library(shiny)
#' ui <- fluidPage(
#'   h3("Validated Numeric Input Examples"),
#'   fluidRow(
#'     column(6, numericInputModuleUI("n1", integerOnly = TRUE, minVal = 0, maxVal = 100, defaultVal = 10)),
#'     column(6, numericInputModuleUI("n2", integerOnly = FALSE, defaultVal = 3.14))
#'   ),
#'   fluidRow(
#'     column(6, numericInputModuleUI("n3", integerOnly = TRUE, minVal = 10)),
#'     column(6, numericInputModuleUI("n4", integerOnly = TRUE, maxVal = 50))
#'   ),
#'   hr(),
#'   verbatimTextOutput("values")
#' )
#'
#' server <- function(input, output, session) {
#'   v1 <- numericInputModuleServer("n1")
#'   v2 <- numericInputModuleServer("n2")
#'   v3 <- numericInputModuleServer("n3")
#'   v4 <- numericInputModuleServer("n4")
#'
#'   output$values <- renderPrint({
#'     list(n1 = v1(), n2 = v2(), n3 = v3(), n4 = v4())
#'   })
#' }
#'
#' shinyApp(ui, server)
#' }
#' @export
numericInputModuleUI <- function(id,
                                 integerOnly = TRUE,
                                 label = NULL,
                                 minVal = NULL,
                                 maxVal = NULL,
                                 defaultVal = NA) {
  ns <- shiny::NS(id)

  baseText <- if (integerOnly) "Enter an integer" else "Enter a value"

  if (!is.null(label)) {
    label_final <- label
  } else {
    label_final <- baseText
    if (!is.null(minVal) && !is.null(maxVal)) {
      label_final <- paste0(label_final, " (", minVal, " ≤ value ≤ ", maxVal, ")")
    } else if (!is.null(minVal)) {
      label_final <- paste0(label_final, " (≥ ", minVal, ")")
    } else if (!is.null(maxVal)) {
      label_final <- paste0(label_final, " (≤ ", maxVal, ")")
    }
  }

  # Only include min/max if not NULL
  args <- list(
    inputId = ns("num"),
    label   = label_final,
    value   = defaultVal
  )
  if (!is.null(minVal)) args$min <- minVal
  if (!is.null(maxVal)) args$max <- maxVal

  tagList(
    do.call(shiny::numericInput, args),
    # hidden input to store integerOnly flag
    shiny::tags$input(
      id = ns("integerOnly"),
      type = "hidden",
      value = if (integerOnly) "1" else "0"
    ),
    shiny::tags$script(HTML(sprintf("
      (function() {
        var el = document.getElementById('%s-num');
        if (!el) return;
        function sendCommitted() {
          var integerOnlyVal = document.getElementById('%s-integerOnly').value;
          Shiny.setInputValue('%s-num_committed',
            { value: el.value, nonce: Date.now(), min: el.min, max: el.max, integerOnly: integerOnlyVal },
            { priority: 'event' });
        }
        el.addEventListener('keydown', function(e) { if (e.key === 'Enter') sendCommitted(); });
        el.addEventListener('blur', sendCommitted);
      })();
    ", id, id, id)))
  )
}


#' Validated Numeric Input Module - Server
#'
#' @param id Module ID.
#' @param invalidMsg Optional custom invalid input message. If provided, it overrides auto-generated messages.
#' @param allowNA Logical. If TRUE (default), blank or NA input is accepted. If FALSE, blank or NA triggers validation error.
#' @param labelMsg Optional string to use in the invalid popup message. If NULL, a message is generated automatically.
#' @return Reactive expression returning the current valid value (`NA` if blank and allowed).
#' @export
numericInputModuleServer <- function(id,
                                     invalidMsg = NULL,
                                     allowNA = TRUE,
                                     labelMsg = NULL) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      prev_value <- shiny::reactiveVal(NA_real_)

      # initialize with the UI's starting value once available
      shiny::observeEvent(input$num, {
        prev_value(input$num)
      }, once = TRUE)

      shiny::observeEvent(input$num_committed, {
        new_str <- input$num_committed$value
        new_val <- suppressWarnings(as.numeric(new_str))

        uiMin <- if (!is.null(input$num_committed$min) && input$num_committed$min != "") as.numeric(input$num_committed$min) else NULL
        uiMax <- if (!is.null(input$num_committed$max) && input$num_committed$max != "") as.numeric(input$num_committed$max) else NULL
        integerOnly <- if (!is.null(input$num_committed$integerOnly) && input$num_committed$integerOnly == "1") TRUE else FALSE

        isNAinput <- is.null(new_str) || new_str == "" || toupper(new_str) == "NA"
        if (isNAinput) {
          if (allowNA) {
            prev_value(NA_real_)
            shiny::updateNumericInput(session, "num", value = NA)
            return()
          } else {
            new_val <- NA_real_
          }
        }

        is_integer <- !is.na(new_val) && (!integerOnly || abs(new_val - round(new_val)) < 1e-9)
        meets_min <- is.null(uiMin) || (!is.na(new_val) && new_val >= uiMin)
        meets_max <- is.null(uiMax) || (!is.na(new_val) && new_val <= uiMax)

        if (!is.na(new_val) && is_integer && meets_min && meets_max) {
          val <- if (integerOnly) as.integer(round(new_val)) else new_val
          prev_value(val)
          shiny::updateNumericInput(session, "num", value = val)
        } else {
          msg <- if (!is.null(invalidMsg)) {
            invalidMsg
          } else if (!is.null(labelMsg)) {
            labelMsg
          } else {
            baseText <- if (integerOnly) "an integer" else "a number"
            tmp <- paste0("Value must be ",
                          if (allowNA) "blank (NA), 'NA', or " else "",
                          baseText)
            if (!is.null(uiMin) && !is.null(uiMax)) {
              tmp <- paste0(tmp, " between ", uiMin, " and ", uiMax)
            } else if (!is.null(uiMin)) {
              tmp <- paste0(tmp, " ≥ ", uiMin)
            } else if (!is.null(uiMax)) {
              tmp <- paste0(tmp, " ≤ ", uiMax)
            }
            tmp
          }
          shiny::showModal(
            shiny::modalDialog(
              title = NULL,
              msg,
              size = "s",
              easyClose = TRUE,
              footer = NULL
            )
          )
          shiny::updateNumericInput(session, "num", value = prev_value())
        }
      }, ignoreInit = TRUE)

      return(shiny::reactive(prev_value()))
    }
  )
}

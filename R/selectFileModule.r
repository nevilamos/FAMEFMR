#' Module UI function for file selection
#' @details Select File Module UI used with selectFileServer
#' modified from https://stackoverflow.com/questions/38747129/how-to-use-shinyfiles-package-within-shiny-modules-namespace-issue
#'
#' @param id string namespace id for module (same as in server)
#' @param label text to appear on select file button eg "Select file"
#' @param title text title for button (not used) default ""
#' @param multiple TRUE for allowing mutiple selection, default FALSE to select single file
#'
#' @return datapath path of selected file ( relative to root_dirs argument in fileSelectServer)
#' @export
#'
#' @examples
#'#simple file selection shiny app
#'library(shiny)
#'ui <-   fluidPage(
#'  selectFileUI(
#'    id = "fileSelect1",label ="Select File"
#'  ),
#'  textOutput("datapath")
#'
#')
#'
#'server <- function(input, output, session) {
#'  rv<-reactiveValues()
#'  x <- selectFileServer(id = "fileSelect1")
#'  observe(rv$datapath<-x$datapath)
#'  observe(output$datapath<-renderText(rv$datapath))
#'
#'
#'}
#'shinyApp(ui,server)

selectFileUI <- function(id, label = "Select file", title ="", multiple = FALSE) {
  # Create a namespace function using the provided id
  ns <- shiny::NS(id)

  shiny::tagList(
    shinyFiles::shinyFilesButton(
      id = ns("get_file_path"),
      label =label,
      title = title,
      multiple = FALSE)
  )

}

#' Module Server function for file selection
#' @details Select File Module UI used with selectFileUI
#' modified from https://stackoverflow.com/questions/38747129/how-to-use-shinyfiles-package-within-shiny-modules-namespace-issue
#'
#' @param id string namespace id for module (same as in server)
#' @param root_dirs named list of paths that for the root directories within
#'  which files may be selected default is the cwd ie c(root = "."))
#' @param filetypes vector of filetypes (by file extension) that may be selected from
#' eg c("txt", "csv", "gpkg")
#'
#' @return filepath as reactiveValue named same as module id  ( relative to root_dirs argument in fileSelectServer)
#' @export
#'
#' @examples
#'#simple file selection shiny app
#'library(shiny)
#'ui <-   fluidPage(
#'  selectFileUI(
#'    id = "fileSelect1",label ="Select File"
#'  ),
#'  textOutput("datapath")
#'
#')
#'
#'server <- function(input, output, session) {
#'  rv<-reactiveValues()
#'  x <- selectFileServer(id = "fileSelect1")
#'  observe(rv$datapath<-x$datapath)
#'  observe(output$datapath<-renderText(rv$datapath))
#'
#'
#'}
#'shinyApp(ui,server)
selectFileServer <- function(id,
                             root_dirs = c(root = "."),
                             filetypes= NULL){
  shiny::moduleServer(id,
               function(input,
                        output,
                        session
               ){
                 rv <- shiny::reactiveValues()
                 shinyFiles::shinyFileChoose(input,
                                             id = "get_file_path",
                                             roots = root_dirs,
                                             session = session,
                                             filetypes = filetypes)
                 datapathtab <-
                   shiny::reactive(
                     shinyFiles::parseFilePaths(
                       roots = root_dirs,
                       input$get_file_path)
                   )
                 shiny::observe(rv$datapath<-datapathtab()$datapath)
                 return(rv)
               }
  )
}

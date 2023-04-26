#' uploadFile UI part of uploadFle shiny module
#'
#' @param id module id name
#' @param label label for shiny fileInput button
#' @param multiple default TRUE allows upload of multiple
#'  files depending on filetype argument in server
#'
#' @return dataframe of files fed to server fileCopy
#' @export
#'
#' @examples

#' #app to upload files into cwd/temp
#'
#' if(!dir.exists("temp")){dir.create("temp")}
#'
#' ui <- fluidPage(
#' uploadFileUI(id = "myID",
#' label = "Select files",
#' multiple = T)
#' )
#'
#' server <- function(input, output) {
#'
#' uploadFileServer(id = "myID",
#' saveToPath = "./temp",
#' filetype ="geog")
#'
#'}
#'shiny::shinyApp(ui = ui, server = server)
uploadFileUI <-  function(id,
                          label = "Select file",
                          multiple = TRUE) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fileInput(
      inputId = ns("files"),
      label = label,
      multiple = multiple,
      options(shiny.maxRequestSize=5*1024^3)
    ),
    shiny::htmlOutput(outputId = ns("myText"))
  )
}



#' uploadFile Server part of uploadFle shiny module
#'
#' @param id module id name
#' @param filetype one of three alternatives:
#' "all" which allows any selected file(s) to be uploaded
#' "geog" which allows load of either a geopackage (.gpkg) file OR
#'  the four critical components of a shapefile ( .shp, .prj, .shx,.dbf) OR
#'  a zip file containing one or more shapefiles, geopackages or
#'  esri geodatabases. Thee zip file is unzipped after upload. OR
#'  any other specified file suffix - tests for and allows upload
#'  only of the defined suffix
#' @param saveToPath Path to save files to on server
#'
#' @return text message on status of upload for printing to ui
#' @export
#'
#' @examples
#' #app to upload files into cwd/temp
#'
#' if(!dir.exists("temp")){dir.create("temp")}
#'
#' ui <- fluidPage(
#' uploadFileUI(id = "myID",
#' label = "Select files",
#' multiple = T)
#' )
#'
#' server <- function(input, output) {
#'
#' uploadFileServer(id = "myID",
#' saveToPath = "./temp",
#' filetype ="geog")
#'
#'}
#'shiny::shinyApp(ui = ui, server = server)
uploadFileServer <-
  function(id,
           filetype = c(NULL,"geog","any other file suffix can be used ( without the dot)"),
           saveToPath = "path to save files to on server")
  {shiny::moduleServer(
    id,
    function(input, output, session) {
      shiny::observeEvent(
        input$files,
        {myInput = input$files
        shapefile_components <-c("shp", "shx", "prj", "dbf")
        L <- length(myInput$name)
        myInput$saveToPaths <-file.path(saveToPath, myInput$name)
        if (L == 0) {myText = ""}
        else if (filetype == "geog") {
          if (L == 1 & all(tools::file_ext(myInput$name) == ("gpkg")))
            {
            file.copy(myInput$datapath,file.path(myInput$saveToPaths),overwrite = T)
            myText <- "geopackage uploaded"
            }
          else if (L == 1 & all(tools::file_ext(myInput$name) == "zip")) {
              file.copy(myInput$datapath,file.path(myInput$saveToPaths),overwrite = T)
            zippedFiles <-unzip(myInput$saveToPaths, list = T)
            correctFiles <-tools::file_ext(
              unlist(
                purrr::map(
                strsplit(
                  zippedFiles$Name,
                  "/"),
                1))) %in% c("gdb", "gpkg", shapefile_components)
            if (all(correctFiles)) {
              unzip(myInput$saveToPaths,exdir = saveToPath,overwrite = T)
              unlink(myInput$saveToPaths)
              myText <-"zipped shapefile(s) and/or gdb or gpkg uploaded and unzipped"
            }else {
              myText <-"<span style=\"color:red\">ERROR <br>
              Incorrect File Selection:<br>
              \nyou have not selected a zip file containing only a <br>
              \ngeopackage ( .gpkg) file or <br>\n
              ESRI geodatabase (.gdb)   or <br>\n
              or one or more components of a shapefile <br>\n
              (.shp,.shx,.dbf,.prj) are missing</span>"

              unlink(myInput$saveToPaths)
            }

          } else if (L == 4 &
                     all(shapefile_components %in%
                         tools::file_ext(myInput$name)) &
                     length(unique(
                       tools::file_path_sans_ext(myInput$name))) ==1) {
            file.copy(myInput$datapath,
                      file.path(myInput$saveToPaths),
                      overwrite = T)
            myText <- "shapefile uploaded"

          } else {
            myText <-"<span style=\"color:red\">ERROR \n
                                you have not selected either a:<br>\n
                                single geopackage ( .gpkg) or a<br>\n
                                zip file containing .gdb .gpkg or shapefiles<br>\n
                                or one or more of .shp,.shx,.dbf,.prj are missing<br>\n
                                or additional incorrect files selected</span>"
          }
        }
        else if (filetype == "all") {
          file.copy(myInput$datapath,
                    file.path(myInput$saveToPaths),
                    overwrite = T)
          myText <- "file(s) uploaded"

        }
        else {
          if (all(tools::file_ext(myInput$name) == (filetype))) {
            file.copy(myInput$datapath,
                      file.path(myInput$saveToPaths),
                      overwrite = T)
            myText <- paste (filetype,"file(s) uploaded")
          } else {
            myText <-
              paste0("<span style=\"color:red\">ERROR \n you
                have selected files with suffix other than'",
                filetype,
                "'</span>")
            }
          }

        output$myText <- shiny::renderText(myText)

        })
            }
    )
      }

#' Uploader for shapefiles zipped geodatabase or geopackage files to shiny app
#'
#' @param myInput  data frame of uploaded files selected using shiny inputFiles()
#' these files may be either:
#'                 a single zip file containing shapefile components (.shp
#'                 .shx,.prj,.dbf) and or one or more esri geodatabases (.gdb)
#'                  and or geopackage files
#'                  or:
#'                  a single geopackage file (.gpkg)
#'                  four files being  all four base components of a shapefile (.shp
#'                 .shx,.prj,.dbfp)
#' @param saveToPath path to save files to on server
#'
#' @return text message indicating success or otherwise of upload for rendering in shiny ui
#' @export
uploadShapeFileGdbGpkg <-
  function(myInput = "data frame of input files from shiny inputFiles",
           saveToPath = "path to save files to on server")
  {
    shapefile_components <- c("shp", "shx", "prj", "dbf")
    L <- length(myInput$name)
    #print(paste("L=",L))
    myInput$saveToPaths <- file.path(saveToPath, myInput$name)
    if (L==0) {
      #print("L==0")
      myText = ""
    } else if (L == 1 & all(tools::file_ext(myInput$name) == ("gpkg"))) {
      file.copy(myInput$datapath,
                file.path(myInput$saveToPaths),
                overwrite = T)
      myText <- "geopackage uploaded"

    } else if (L == 1 & all(tools::file_ext(myInput$name) =="zip"))  {
      file.copy(myInput$datapath,
                file.path(myInput$saveToPaths),
                ,overwrite = T)
      ###insert check for gdb in zip if correct then unzip it, if not delete file and return error
      zippedFiles <- unzip(myInput$saveToPaths, list = T)
      #print("zippedFiles")
      correctFiles <-
        tools::file_ext(unlist(purrr::map(strsplit(
          zippedFiles$Name, "/"
        ), 1))) %in% c("gdb", "gpkg", shapefile_components)
      if (all(correctFiles)) {
        #print("all correct")
        unzip(myInput$saveToPaths,exdir = saveToPath,overwrite = T)
        #print("unzipped file")
        unlink(myInput$saveToPaths)
        myText <-
          paste("zipped shapefile(s) and/or gdb or gpkg uploaded and unzipped")
      } else if (sum(correctFiles) > 0) {
        #print("some correct")
        unzip(myInput$saveToPaths, files = zippedFiles$Name[correctFiles],exdir = saveToPath,overwrite = T)
        #print("unzipped file")
        unlink(myInput$saveToPaths)
        myText <-
          paste("zipped shapefile(s) and/or gdb or gpkg uploaded and unzipped")
      } else{
        #print("bailed out")
        myText <-
          paste(
            "<span style=\"color:red\">ERROR <br> Incorrect File Selection:<br>
            you have not selected a zip file containing a <br>
            geopackage ( .gpkg) file or <br>
            ESRI geodatabase (.gdb) file  or <br>
            or one or more componets of a shapefile <br>
            (.shp,.shx,.dbf,.prj) are missing</span>"
          )
        unlink(myInput$saveToPaths)

      }


    } else if (L == 4 &
               all(shapefile_components %in% tools::file_ext(myInput$name)) &
               length(unique(tools::file_path_sans_ext(myInput$name))) == 1) {
      file.copy(myInput$datapath,
                file.path(myInput$saveToPaths),
                overwrite = T)
      myText <- "shapefile uploaded"
    } else {
      myText <-
        "<span style=\"color:red\">ERROR \n your have not selected either a:<br>
      single geopackage ( .gpkg) or a<br>
      zip file containing .gdb .gpkg or shapefiles<br>
      or one or more of .shp,.shx,.dbf,.prj are missing<br>
      or additional incorrect files selected</span>"
    }

    myText
  }

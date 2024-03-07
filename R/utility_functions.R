
############################################################################
#Utility functions used in FAME Fire History analyses (often called within the
#main analysis functions) written by nevil.amos@delwp.vic.gov.au
############################################################################


#' remove empty directories from path
#' @details Removes all empty subdirectories from the nominated path does not remove the nominated path directory even if empty
#' @param rootDir relative path to remove all empty subdirectories from. Default
#'   value "./Results"
#' @export
removeEmptyDirs <- function(rootDir="./Results"){

  dirList <- list.dirs(rootDir)[-1] #[-1] makes sure the root results directory is not deleted
  if(length(dirList) > 0){
    delList <- dirList[sapply(dirList, function(x) length(list.files(x)) == 0)]
    while((length(delList) > 0) & (length(dirList) > 0)){
      unlink(delList, recursive = TRUE)
      dirList <- list.dirs(rootDir)[-1]
      delList <- dirList[sapply(dirList, function(x) length(list.files(x)) == 0)]
    }
  }
}


#' Adds concatenated String of X and Y coordinates of centroids of polygons to
#' Simple Features polygon object. This String acts as a key to identify
#' spatially identical polygons for use in tidyverse pivot functions.
#' @param myDF sf polygon object
#' @return character vector of XYStrings
#' @export
add_xystring <- function(myDF){
  Coords <- as.data.frame(
    do.call(
      rbind,
      sf::st_geometry(sf::st_centroid(myDF))))
  names(Coords) <- c("X","Y")
  XYString <- paste(Coords$X, Coords$Y, sep = "")
  x <- dplyr::mutate(myDF, XYString)
  return(x)
}

#' Extract VBA (Victorian Biodiversity Atlas)  species ID numbers from file
#' paths extracts four or five digit species numbers (Victorian Biodiversity
#' Atlas TAXON_IDs) from vector of paths or file names containing files of e.g.
#' species HDMS containing the 5 digit TAXON_ID in their name
#' @param x Vector of species file Pathnames containing VBA numbers
#' @return numeric vector of 4or 5 digits (ususally TAXON_ID)
#' @export
get_Spp_No <- function( x = "Vector of Sp file Pathnames"){
  Fnames= basename(x)
  pos = regexpr("[0-9][0-9][0-9][0-9]+", Fnames)
  myEnd = pos-1+attr(pos, "match.length")
  y = as.numeric(substr(Fnames, pos, myEnd))
  return(y)
}

#' Joins one or more lookup tables to table containing ID values Function joins
#' Lookup tables (LUTS) to dataframe containing ID_NO: Name combinations
#'
#' @param myDF dataframe or similar containing indices for the LUTS listed, to
#'   which the LUTS will be dplyr::left_joined
#' @param LUTS vector of names of LUTS in memory defaults
#'   =c("TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")
#' @return a data.frame with the LUTS joined to it
#' @export
Join_Names <- function(myDF,   #dataframe or similar containing indices for the LUTS listed
                       LUTS = c("TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")){
  for(i in LUTS){
    try(myDF <- dplyr::left_join(myDF,get(i)))}
  return (myDF)
}


#' Calculates multiplier to convert from raster cell count to area in hectares
#' @param RasterMetres numeric Value cell resolution in Metres (usually from
#'   RasterRes in settings file).
#' @return numeric Multiplier to convert cell count to area in hectares
#' @export
cellsToHectares <- function(RasterMetres = RasterRes){
  (RasterMetres / 100) ^ 2
}


#' Checks whether all values in  one vector are in another vector
#' @param x Vector of values to check if all are in second vector
#' @param v Second vector of values that may or may not contain all values in x
#' @return logical
#' @export
notAllIn <- function(x, v = V){
  anyNA(unlist(lapply(x, match, v)))
}





#' Fix Pivot_wider list of lists columns
#' @details Supporting function to deal with pivot_wider returning list of lists in some cases
#' @param df wide format data frame with fields that are lists of lists
#'
#' @return wide format data frame without fields that are lists of lists
#' @export

unlistPivot_wider <- function(df){
  df[unlist(lapply(df , is.null))] <- NA
  Y <- unlist(df)
  return (Y)
}



#' Save settings to csv file
#'
#' untility helper function for FAME to make easily readable file of settings used
#' in analysis,stored with results output.
#'
#' @param outputSettings list of settings accumaulated in FAME processing
#' @param resultsDir the resuts directory where the ouptuSettings.csv file is saved
#'
#' @return writes a table as csv file to resultsDir
#' @export
saveSettingsTable<-function(outputSettings,resultsDir){
  # convert list to dataframe with each item as row
  outTab<-as.data.frame(t(as.data.frame((outputSettings))))
  # name columns
  names(outTab)[1]<-"Value"
  write.csv(outTab,file.path(resultsDir,"outputSettings.csv"))

}

#' Save Settings for fhProcess()
#'
#' utility/helper function to capture the settings used in creating the fire
#' history (fhProcess()) analysis part of FAME
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#'
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#' @export
#'
saveFHsettings<-function(outputSettings,rv){
  outputSettings$rawFHPath[1] = rv$rawFHPath
  outputSettings$resultsDir[1] = rv$resultsDir
  outputSettings$outputFH[1] = rv$outputFH
  outputSettings$RasterRes[1] = rv$FHAnalysis$RasterRes
  outputSettings$ClipPolygonFile[1] =  rv$FHAnalysis$ClipPolygonFile
  outputSettings$REGION[1] = rv$FHAnalysis$Region_No
  outputSettings$public[1] = rv$FHAnalysis$PUBLIC_ONLY
  outputSettings$Start_Season[1] = rv$startTimespan
  outputSettings$fhAnalysisName [1] =rv$FHAnalysis$name
  return(outputSettings)
}

#' Save Settings for PU file
#'
#' utility/helper function to capture the settings used in FAME analysis
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#' @export
#'
savePUsettings<-function(outputSettings,rv){
  #outputSettings$usePUPolys = input$usePUpolys
  outputSettings$puPath = rv$puPath
  return(outputSettings)
}

#' Save Settings for species calcuations
#'
#' utility/helper function to capture the settings used in FAME analysis
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#'
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#' @export
#'
saveSPsettings<-function(outputSettings,rv){
  outputSettings$fhAnalysisName [1] =rv$FHAnalysis$name
  outputSettings$startBaseline[1] = rv$startBaseline
  outputSettings$endBaseline[1] = rv$endBaseline
  outputSettings$spListChoice[1] = rv$spListChoice
  outputSettings$TaxonListPath[1] = rv$TaxonListPath
  outputSettings$spResponseChoice [1] =rv$spResponseChoice
  outputSettings$SpGSResponses[1] = rv$SpGSResponses
  outputSettings$WriteSpRasters[1] = rv$makeRArasters
  outputSettings$yearsForRasters[1] = paste(rv$yearsForRasters,collapse=' ')
  return(outputSettings)
}



#' Save JFMP settings
#'
#' utility/helper function to capture the settings used in FAME analysis
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#'
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#' @export
saveJFMPsettings<-function(outputSettings,rv){
 outputSettings$ranJfmp[1] = TRUE
 outputSettings$myJFMPSeason0[1] = rv$JFMPSeason0
 outputSettings$puPath[1] = rv$puPath
 outputSettings$zoneWtFile[1] = rv$zoneWtFile
 outputSettings$jfmpMetricWtFilePath[1] = rv$jfmpMetricWtFilePath
 outputSettings$targetHaFilepath[1] = rv$targetHaFilepath
 return(outputSettings)
}


#' Save settings to csv file
#'
#' utility helper function for FAME to make easily readable file of settings used
#' in analysis,stored with results output.
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param outputNames data.frame formatted and named as the settingstable will be output, but with column 2 values as NA
#' @param resultsDir the results directory where the ouptutSettings.csv file is saved
#'
#' @return writes a table as csv file to resultsDir
#' @export
#'
#'
saveSettingsTable<-function(outputSettings,resultsDir,outputNames){
  # convert list to dataframe with each item as row
  outTab<-as.data.frame(t(as.data.frame((outputSettings))))
  names(outTab)<-names(outputNames)
  write.csv(outTab,file.path(resultsDir,"outputSettings.csv"),row.names = TRUE,
            col.names = TRUE)

}

#' Save Settings for fhProcess()
#'
#' utility/helper function to capture the settings used in creating the fire
#' history (fhProcess()) analysis part of FAME
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#'
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#' @export
#'
saveFHsettings<-function(outputSettings,rv){
  outputSettings$rawFHPath[1] = rv$rawFHPath
  outputSettings$resultsDir[1] = rv$resultsDir
  outputSettings$outputFH[1] = rv$outputFH
  outputSettings$RasterRes[1] = rv$FHAnalysis$RasterRes
  outputSettings$max_interval[1] = rv$FHAnalysis$max_interval
  outputSettings$ClipPolygonFile[1] =  rv$FHAnalysis$ClipPolygonFile
  outputSettings$REGION[1] = rv$FHAnalysis$Region_No
  outputSettings$public[1] = rv$FHAnalysis$PUBLIC_ONLY
  outputSettings$Start_Season[1] = rv$startTimespan
  outputSettings$fhAnalysisName [1] =rv$FHAnalysis$name
  return(outputSettings)
}

#' Save Settings for PU file
#'
#' utility/helper function to capture the settings used in FAME analysis
#'
#' @param outputSettings list of settings accumulated in FAME processing
#'
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#' @export
#'
savePUsettings<-function(outputSettings,rv){
  #outputSettings$usePUPolys = input$usePUpolys
  outputSettings$puPath = rv$puPath
  return(outputSettings)
}

#' Save Settings for species calcuations
#'
#' utility/helper function to capture the settings used in FAME analysis
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#'
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#' @export
#'
#'
saveSPsettings<-function(outputSettings,rv){
  outputSettings$fhAnalysisName [1] =rv$FHAnalysis$name
  outputSettings$startBaseline[1] = rv$startBaseline
  outputSettings$endBaseline[1] = rv$endBaseline
  outputSettings$spListChoice[1] = rv$spListChoice
  outputSettings$TaxonListPath[1] = rv$TaxonListPath
  outputSettings$spResponseChoice [1] =rv$spResponseChoice
  outputSettings$SpGSResponses[1] = rv$SpGSResponses
  outputSettings$WriteSpRasters[1] = rv$makeRArasters
  outputSettings$yearsForRasters[1] = paste(rv$yearsForRasters,collapse=' ')
  return(outputSettings)
}



#' Save JFMP settings
#'
#' utility/helper function to capture the settings used in FAME analysis
#'
#' @param outputSettings list of settings accumulated in FAME processing
#' @param rv reactivevalues() object including settings defined in shiny app or
#' list of same object made in script version
#'
#' @return outputSettings list of settings accumulated in FAME processing with
#' updated values
#'
#' @export
#'
saveJFMPsettings<-function(outputSettings,rv){
  outputSettings$ranJfmp[1] = TRUE
  outputSettings$myJFMPSeason0[1] = rv$JFMPSeason0
  outputSettings$puPath[1] = rv$puPath
  outputSettings$zoneWtFile[1] = rv$zoneWtFile
  outputSettings$jfmpMetricWtFilePath[1] = rv$jfmpMetricWtFilePath
  outputSettings$targetHaFilepath[1] = rv$targetHaFilepath
  return(outputSettings)
}

#
#' Align raster extent with template grid
#'
#' @param inputPath the path to the raster to be processed to match template grid
#' @param RES the Resolution of the template grid ( units depend on CRS) for FAME it will be 75 0r 225 metres
#' @param CRS in form acceptable for terra::crs() function of   default VicGrid94 given as "epsg:3111"
#' @param EXTENT extent of the template raster as vector giving xmin ,ymin xmax,ymax in same units as resolution.
#' @param outPath path to ouptut file to (folder path) if null then file is output to same directory as input with RES_cutom.tif appended to filename
#'
#' @return outputRaster of same extent, resolution and CRS as the template settings
#' @export
#'
#' @examples
#'  terra::writeRaster(terra::rast(matrix(10,10,10),ext=c(143,143.3,-36,-35.7),crs="epsg:4283" ),"temp.tif",overwrite=TRUE)
#'  alginCustomRaster(inputPath ="./temp.tif", RES =  75, CRS="epsg:3111",EXTENT = c(2126550, 2941050, 2258825, 2828075 ))
#'  unlink("./temp_75_custom.tif")
#'  unlink("temp.tif")
alginCustomRaster<- function(inputPath ="./MEWCustom/MEW_new.tif", RES =  75, CRS="epsg:3111",EXTENT = c(2126550, 2941050, 2258825, 2828075 ),outPath = NULL){


  template = terra::rast(crs=CRS,extent = EXTENT,res = RES)
  inputRaster<-terra::rast(inputPath)
  print(" INPUT:")
  print( inputRaster)
  if (terra::crs(inputRaster)==""){
    stop("input raster does not have a coordinate reference system/ projection")
  }else if (identical(terra::crs(inputRaster),terra::crs(template))){
    outputRaster<- terra::resample(inputRaster,template,method = "near")
    print("resampled")
  }else {
    outputRaster<-terra::project(inputRaster,template,method = "near")
    print("reprojected")
  }
  if (is.null (outPath)){
    outPath<-paste0(tools::file_path_sans_ext(inputPath),"_",RES,"_custom.tif")}
  else{outPath<-outPath}

  #raster is written out to the same directory as the input with the resolution pasted into its name
  terra::writeRaster(outputRaster,outPath,overwrite=T)
  print(" OUTPUT:")
  print(outPath)
  return(outputRaster)
}




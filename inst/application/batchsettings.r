#Settings file to run FAME analysis from script without shinyApp
#path(s) to your rawFH file ( output of the ARCGIS preprocessing tool) each given as a separate string of the full path or ./path relative to the working directory. 

#the root directory in which results directories should be created and outputs saved

resultsDir = "./results"

#fhProcess settings
#path to the input fire history geodatabase (".gdb" or ".gpkg") or shapefile (".shp") usually provided in settings
firstFH	="rawFH/FH.gpkg"

#if rawFH is a geodatabase the name of the layer containing the Fire History if this is not provided and the raw FH is a geodatabase or geopackage then the first layer in the file is returned.This argument is available only in the r script version - it is not exposed in the shinyApp where it is set to NULL

firstFHLayer	=NULL
#Second fire history to be combined with FH1 to make a fire scenario same formats as for firstFH or NULL if there is only a single FH used.
secondFH	=NULL

#if secondFH is a geodatabase the name of the layer containing the Fire History if this is not provided and the raw FH is a geodatabase or geopackage then the first layer in the file is returned.This argument is available only in the r script version - it is not exposed in the shinyApp where it is set to NULL

secondFHLayer	=NULL

#value to use for cases where fire type is "OTHER" or "UNKNOWN",1 ="BURN",2="BUSHFIRE",NA = Fire excluded from analysis default is 2 ("BUSHFIRE")
otherUnknown = 2

#Default NULL otherwise four digit integer SEASON for fire applied across the whole bounding box
baseFire	=1755

 #integer number of years maximum inter fire interval after a HIGH fire for which subsequent fire is reassigned to HIGH default is NULL in which case there is no reassignment
max_interval	=  5#


startTimespan	=NA

#integer Last SEASON required, if NA then largest value in fire history scenario used
endTimespan =	NA

#vector of valid names in the input FIRETYPE column in the input fire history dataset(s), if the column contains NA or values not on this list an error will occur
validFIRETYPE	 = c("BURN", "BUSHFIRE", "OTHER", "UNKNOWN")


#the integer value of the region number 1-6 for FFR regions,7 for user supplied adHoc polygon, 99 for Statewide
REGION_NO = "7"  #"99"

#path to the ad Hoc polygon if REGION_NO 			== 7 otherwise NULL
AdHocPath = "AdHocPolygons/AA.shp"

doSpeciesCalculations= TRUE

# whether or not custom species list is used ( TRUE for custom species list) otherwise FALSE
spListChoice = TRUE
customSpList = "CustomCSV/5spWithAreas.csv"

#whether the species response file contains abundance data by 4 growth stages ( TRUE) or YSF (FALSE) 0-400
abundByGS =TRUE
# whether the default or custom response file is used for abundance responses
spResponseChoice = FALSE
#path of custom response abundance file if spResponseChoice==TRUE
customResponseFile = #"C:/Data/FAMEshiny/jennyOtawysSbbdebug/sbb_gam1_fame_test.csv"


#whether analysis should be undertaken only on public land ( public land only if TRUE)
public = TRUE
#raster resolution 75 or 225
RasterRes = 225

#"first season for which output is wanted ( four digit year as integer)
#startTimespan = 1988
#endTimespan =2030
#start and end baseline years if single year then the two values should be equal
startBaseline = 2000
endBaseline = 2010
#endSEASON = NULL


#whether to write TFI rasters
makeTFIrasters = TRUE
#whether to write BBTFI rasters
makeBBTFIrasters = TRUE

#whether to write Species Relative Abundance rasters
makeRArasters = TRUE


# which years to write rasters if makeRArasters = TRUE
yearsForRasters = NULL#  c(2020)vector of SEASONS within limit of , or NULL for all years


#whether a PUpolys file is to be used - this is required for autoJFMP analysis, and can also be used as a work around to allow grouping of outputs e.g. by fire district for the FAME dashboard ( by providing a PUPolys file where the PU correspond to fire districts rather than  burn units "LF_DISTRICT_with_PU_field".)
usePUpolys = FALSE
puPath = "ReferenceShapefiles/LF_DISTRICT_with_PU_field.shp"
#Settings if JFMP analysis is also to be run this is not working yet so  keep runJFMP = FALSE and JFMPSeason0 = endTimespan
runJFMP = FALSE


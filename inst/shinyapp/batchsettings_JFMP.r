#Settings file to run a number of alternative fire history scenarios for the
#same area of interest though fame as a batch.

# the batch will run through each of the raw fire history scenarios specified
# by the shapefile paths in rawFHPath  in series putting the results in a
# separate subdirectory of the resultsDir.  All other inputs are constant
# throughout the batch

# sets default directory for FAME input files
options(FAME.data_dir = "C:/Data/FAMEshiny-4.6")

#path(s) to your rawFH file ( output of the ARCGIS preprocessing tool) each
#given as a separate string of the full path or ./path relative to the working
#directory.


rawFHPaths = c(fame_data_path("rawFH","HUME2020JFMP_FH_RAW.shp"))
#the root directory in which results directories should be created and outputs saved
<<<<<<< HEAD
resultsDir = fame_data_path("results")
=======
resultsDir = "./results"
>>>>>>> 79aa0dbf8b8144c1c88960897232186556e9fdab
#the integer value of the region number 1-6 for FFR regions,7 for user supplied adHoc polygon, 99 for Statewide
REGION_NO = "4"  #"99"

#path to the ad Hoc polygon if REGION_NO 			== 7 otherwise NULL
AdHocPath = fame_data_path("AdHocPolygons","DemoAdHocPolygon.shp")

doSpeciesCalculations = TRUE

# whether or not custom species list is used ( TRUE for custom species list) otherwise FALSE
spListChoice = TRUE
customSpList = fame_data_path("CustomCSV", "Demo5SpeciesList.csv")

#whether the species response file contains abundance data by 4 growth stages ( TRUE) or YSF (FALSE) 0-400
abundByGS = TRUE
# whether the default or custom response file is used for abundance responses
spResponseChoice = FALSE
#path of custom response abundance file if spResponseChoice==TRUE
customResponseFile = NULL

  #value to use for cases where fire type is "OTHER" or "UNKNOWN",1 ="BURN",2="BUSHFIRE",NA = Fire excluded from analysis default is 2 ("BUSHFIRE")
  otherAndUnknown = 2
#whether analysis should be undertaken only on public land ( public land only if TRUE)
public = TRUE
#raster resolution 75 or 225
RasterRes = 225

#"first season for which output is wanted ( four digit year as integer)
startTimespan = 1980
#"first season for which output is wanted ( four digit year as integer)
endTimespan = 2024
#start and end baseline years if single year then the two values should be equal
startBaseline = 1980
endBaseline = 1980
endSEASON = NULL


#whether to write TFI rasters
makeTFIrasters = TRUE
#whether to write BBTFI rasters
makeBBTFIrasters = TRUE
#whether to write growth Stage rasters
makeGSrasters =TRUE
#whether to write Species Relative Abundance rasters
makeRArasters = TRUE


# which years to write rasters if makeRArasters = TRUE
yearsForRasters = NULL#  c(2020)vector of SEASONS within limit of , or NULL for all years


#whether a PUpolys file is to be used - this is required for autoJFMP analysis, and can also be used as a work around to allow grouping of outputs e.g. by fire district for the FAME dashboard ( by providing a PUPolys file where the PU correspond to fire districts rather than  burn units "LF_DISTRICT_with_PU_field".)
usePUpolys = TRUE
puPath = "C:/Data/FAMEshiny-4.6/PUPolygons/HUME_PU.shp"

#Settings if JFMP analysis--------------- is also to be run this is not working yet so  keep runJFMP = FALSE and JFMPSeason0 = endTimespan
runJFMP = TRUE
# if JFMP is true all the below settings are required
JFMPSeason0 = 2020
# JFMP Area Target file  ----
targetHaFilepath = "C:/Data/FAMEshiny-4.6/CustomCSV/Hume/hume_JFMPTargets.csv"
# jfmp metric wt File be run ----
jfmpMetricWtFile = "C:/Data/FAMEshiny-4.6/CustomCSV/Hume/hume_JFMPMetricWt.csv"
# zone wt File be run ----
zoneWtFile = "C:/Data/FAMEshiny-4.6/CustomCSV/Hume/hume_JFMPZoneWt.csv"
jfmpSpeciesList = "C:/Data/FAMEshiny-4.6/CustomCSV/Hume/Hume5spwithcellsinarea.csv"





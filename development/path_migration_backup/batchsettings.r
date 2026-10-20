#Settings file to run a number of alternative fire history scenarios for the
#same area of interest though fame as a batch.

# the batch will run through each of the raw fire history scenarios specified
# by the shapefile paths in rawFHPath  in series putting the results in a
# separate subdirectory of the resultsDir.  All other inputs are constant
# throughout the batch

#path(s) to your rawFH file ( output of the ARCGIS preprocessing tool) each
#given as a separate string of the full path or ./path relative to the working
#directory.
rawFHPaths = c("rawFH/rawFH_demoFH2025_vg94_DemoAdHocPolygon_1755.gpkg")
#the root directory in which results directories should be created and outputs saved
resultsDir = "./results"
#the integer value of the region number 1-6 for FFR regions,7 for user supplied adHoc polygon, 99 for Statewide
REGION_NO = "7"  #"99"

#path to the ad Hoc polygon if REGION_NO 			== 7 otherwise NULL
AdHocPath = "AdHocPolygons/DemoIrregularPolygon.gpkg"

doSpeciesCalculations = TRUE

# whether or not custom species list is used ( TRUE for custom species list) otherwise FALSE
spListChoice = TRUE
customSpList = "CustomCSV/Demo5SpeciesList.csv"

#whether the species response file contains abundance data by 4 growth stages ( TRUE) or YSF (FALSE) 0-400
abundByGS = TRUE
# whether the default or custom response file is used for abundance responses
spResponseChoice = FALSE
#path of custom response abundance file if spResponseChoice==TRUE
customResponseFile = #"C:/Data/FAMEshiny/jennyOtawysSbbdebug/sbb_gam1_fame_test.csv"
  
  #value to use for cases where fire type is "OTHER" or "UNKNOWN",1 ="BURN",2="BUSHFIRE",NA = Fire excluded from analysis default is 2 ("BUSHFIRE")
  otherUnknown = 2
#whether analysis should be undertaken only on public land ( public land only if TRUE)
public = TRUE
#raster resolution 75 or 225
RasterRes = 225

#"first season for which output is wanted ( four digit year as integer)
startTimespan = 1980
endTimespan = 2030
#start and end baseline years if single year then the two values should be equal
startBaseline = 2000
endBaseline = 2010
endSEASON = NULL


#whether to write TFI rasters
makeTFIrasters = TRUE
#whether to write BBTFI rasters
makeBBTFIrasters = TRUE
#whether to write growth Stage rasters
makeGSrasters =TRUE
#whether to write Species Relative Abundance rasters
makeRArasters = FALSE


# which years to write rasters if makeRArasters = TRUE
yearsForRasters = NULL#  c(2020)vector of SEASONS within limit of , or NULL for all years


#whether a PUpolys file is to be used - this is required for autoJFMP analysis, and can also be used as a work around to allow grouping of outputs e.g. by fire district for the FAME dashboard ( by providing a PUPolys file where the PU correspond to fire districts rather than  burn units "LF_DISTRICT_with_PU_field".)
usePUpolys = FALSE
puPath = "./ReferenceShapefiles/LF_DISTRICT_with_PU_field.shp"
#Settings if JFMP analysis is also to be run this is not working yet so  keep runJFMP = FALSE and JFMPSeason0 = endTimespan
runJFMP = FALSE

JFMPSeason0 = endTimespan

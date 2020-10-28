## -------------------------------------------------------------------------
## R Default Settings for FAMEv1.0.1 batch process--------------------------
# At the moment this does not run the aspatial GSO tool, it runs all the
# spatial TFI and GSO calcuations. 
# if customSppList and customAbundanceLU = NULL then the default files are used. 
# if REGION_NO = 7 then an adHocPolygon shapefile must also be provided (inputs directory)
## -------------------------------------------------------------------------

## DEFAULT SETTINGS --------------------------------------------------------
## DEFAULT REGION NUMBER ---------------------------------------------------
# see look up table REG_LUT for value integers
# 1-6 FFR regions, 7 adHoc polygon (user input), 99 statewide
REGION_NO = 7
# default path to the adHoc polygon (REGION_NO == 7)
adHocPolygon 		= "./AdHocPolygons/DemoAdHocPolygon.shp"

## RAW FIRE HISTORY FILE ---------------------------------------------------
# path to rawFH file (output of the ARCGIS preprocessing tool)
rawFH = "./rawFH/FRAU_FH_2020_DemoAdHocPolygon.shp"

## CUSTOM SPECIES LIST -----------------------------------------------------
#path to custom species list if used
customSpList = "./CustomCSV/DemoCustomSpeciesList.csv"

## RASTER SETTINGS ---------------------------------------------------------
# raster resolution 75 or 225                                               
RasterRes = 225

# path to HDMVals225 file
#HDMVals225	=	"./HDMS/HDMVals225.rdata"                                     #####-----DELETE-----#####

# whether to output species rasters ("Yes" or "No")                         #####-----should this be no? seems odd to say yes here, and null below-----#####
writeSpRasters = "Yes"

# if writeSpRasters = "Yes", vector of years for which rasters are to be written
yearsForRasters = NULL

# whether to output GS rasters ("Yes" or "No")
writeGSRasters = "Yes"

# whether to write TFI rasters ("Yes" or "No")
makeTFIRasters = "Yes"

# whether to output bbtfi rasters ("True" or "False")                       #####-----make yes or no for consistency-----#####
makeBBTFIrasters = TRUE

## SEASON RELATED SETTINGS -------------------------------------------------
# first season for which output is wanted ( four digit year as integer)
# if NUll, second season in in history is used. cannot use first season
# because it has no interval (this may still fail if there is no overlap)
start.SEASON = 2000                                                         #####-----are there user inputs for this, seems hard to change in here----#####
#end.SEASON = 2000                                                          #####-----there is no default end.SEASON defined.

# start and end baseline years
# for calcuation of change in abundance relative to baseline
# if single year then the two values should be equal
startBaseline = 2000
endBaseline = 2000



## SETUP DIRECTORY PATH FOR ANALYSIS OUTPUTS--------------------------------
ResultsDir <- file.path("./Results",
                        paste(tools::file_path_sans_ext(basename(rawFH)), RasterRes, sep = "_")
                        )
dir.create(ResultsDir)
# copy across default settings for record keeping                            #####-----is this for record keeping? should it be for these parameters at the end if they change? E.G. LINE 5 THIS SCRIPT
file.copy("settings.r",
          file.path(ResultsDir, "settings.r")
          )


## BASIC DEFAULT SETTINGS --------------------------------------------------
# path to custom response list if used
customResponseFile = NULL

# Valid values for rawFH FIRETYPE
# see look up table FIRETYPE_LUT
validFIRETYPE = c("BURN","BUSHFIRE","UNKNOWN","OTHER")

# Firetype value (2,1,NA)
# value to use for cases where fire type is "OTHER" or "UNKNOWN"
# see look up table FIRETYPE_LUT (1 ="BURN",2="BUSHFIRE")
# NA = Fire excluded from analysis default is 2 ("BUSHFIRE")                 #####-----should add NA to lut for consistency? (even if not used)-----#####
OtherAndUnknown = 2  # (2,1,NA)                                              #####-----this appears to include other and unknown firetypes as bushfires-----#####

# whether analysis should be undetaken only on public land
PUBLIC_LAND_ONLY = FALSE                                                     #####-----consistency throughout - True/False, yest/no-----#####

# import 4gs table from csv (contains 4gs info by years since fire)
myEFG_TSF_4GS = read.csv("./ReferenceTables/EFG_TSF_4GScorrectedAllEFGto400yrsV2.csv")[,c('EFG_NO', 'GS4_NO', "YSF")]

# parallel processing information
Ncores = 4
print(paste("Using", Ncores, "cores"))

# change below to TRUE for abrreviated test (on first two years only)
Test = FALSE



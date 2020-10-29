## CLEAR WORKSPACE ---------------------------------------------------------#####-----move to RunWholeFAME.r?-----#####
rm(list=ls(all=TRUE))
gc()

## LOAD LIBRARIES AND CUSTOM FUNCTIONS--------------------------------------
options(stringsAsFactors = FALSE)
library(data.table)
library(doParallel)
library(tidyverse)
library(fasterize)
library(filematrix)
library(foreach)
library(Matrix.utils)
library(raster)
library(Rfast)
library(sf)
library(tabularaster)
library(tcltk)
library(tictoc)
library(tools)

# loads functions used in TFI and RA calculations
source("Utility_Functions.R")# brief description  
source("EcoResFunctionsFMRv2.r")  # brief description                       #####-----comment-----#####
source("TFI_functionsFMRv2.r")  # brief description                         #####-----comment-----#####
source("GS_Calcs.R")  # brief description                                   #####-----comment-----#####
source("calc_U_AllCombs.r")  # brief description                            #####-----comment-----#####
source("calcBBTFI_2.R")
                              #####-----comment-----#####
#source("ButtonDisableHelpers.r")                                           #####-----delete------#####

## SET THE MAXIMUM SIZE OF FILES FOR UPLOAD / DOWNLOAD----------------------#####-----delete------#####
#options(shiny.maxRequestSize=2*1024^3)                                     #####-----delete------#####

## START A TIMER FOR WHOLE ANALYSIS-----------------------------------------
tic("whole process time")

## LOAD INITIAL SETTINGS, INPUT FILES AND TABLES----------------------------
tic("initial setup")
source("./settings.r") # imports initial settings and filepaths

## MAKE RESULTS DIRECTORIES ------------------------------------------------
# create a unique results directory for each run of scenarios 
dir.create(ResultsDir)                                                      #####----- does ResultsDIr ever have more than one dir?
                                                  #####----- ^^^  This is done in settings. Consider deleting if not affecting shiny interface
for (d in c("RA_Rasters", "TFI_Rasters", "GS_Rasters", "BBTFI_Rasters"))
  {dir.create(file.path(ResultsDir, d))}
rm(d)
                                                                            #####-----possibly delete all below. test first
                                                                            #####----- does ResultsDIr ever have more than one dir?
# using starting time as a string for Results directory name                #####-----can't see where the starttime is incorporated
#WD<-getwd() # set working directory                                        #####-----possibly delete (WD not called anywhere)-----#####

#cleans up empty results directoris, and creats new                         #####-----delete------#####
#removeEmptyDirs(rootDir = "./")                                            #####-----delete------#####

# for (d in c(ResultsDir)){dir.create(d)}
# rm(d)
# dir.create(file.path(ResultsDir,"RA_Rasters"))
# dir.create(file.path(ResultsDir,"TFI_Rasters"))
# dir.create(file.path(ResultsDir,"GS_Rasters"))
# dir.create(file.path(ResultsDir,"BBTFI_Rasters"))

## MAKE LOOK UP TABLES-----------------------------------------------------#####-----rename all lut to lut_ (group)
# Fire region, state or adhoc polygon for analysis look up table
REG_LUT <- tibble(FIRE_REG = c(99, 1, 2, 3, 4, 5, 6, 7),
                  FIRE_REGION_NAME = c("WHOLE OF STATE",
                                       "BARWON SOUTH WEST",
                                       "GIPPSLAND",
                                       "GRAMPIANS",
                                       "HUME",
                                       "LODDON MALLEE",
                                       "PORT PHILLIP",
                                       "USER DEFINED POLYGON")
                  )

# Fire FMZ look up table
FIREFMZ_LUT <- tibble(FIREFMZ = c(0, 1, 2, 3, 4, 5),
                      FIRE_FMZ_NAME = c("0 - Non FPA Land",
                                        "1 - Asset Protection Zone",
                                        "2 - Bushfire Moderation Zone",
                                        "3 - Landscape Management Zone",
                                        "4 - Planned Burn Exclusion Zone",
                                        "5 - Unknown. Contact Fire Management Officer"),
                      FIRE_FMZ_SHORT_NAME = c("Non FPA",
                                              "APZ",
                                              "BMZ",
                                              "LMZ",
                                              "PBEZ",
                                              "UNK")
                      )

# DELWP region names look up table                                         #####-----any reason these numbers dont match fire region
DELWP_LUT <- tibble(DELWP = c(2,4,3,5,6,7),                                #####-----LUT and then you could just slice that LUT.
                    DELWP_REGION = c("GIPPSLAND", "HUME", "PORT PHILLIP",
                                     "BARWON SOUTH WEST", "LODDON MALLEE", "GRAMPIANS")
                    )

# TFI_STATUS values lookup table and export to csv
TFI_STATUS_LUT <- structure(list(TFI_VAL = c(-99, 0, 1, 5, 6),
                                 TFI_STATUS = c("NONE","WITHIN_TFI", "BELOW_MIN_TFI",
                                                "ABOVE_MAX_TFI", "ABOVE_MAX_BELOW_MIN_HI_TFI")
                                 ),
                            row.names = c(NA,-5L),
                            class = c("tbl_df", "tbl", "data.frame")
                            )
write.csv(TFI_STATUS_LUT,
          file.path(ResultsDir, "TFI_Rasters", "TFI_STATUS_LUT.csv")
          )

# Firetype look up table
FIRETYPE_LUT <- tibble(TYPE = c(1, 2),
                       FIRETYPE = c("BURN", "BUSHFIRE")
                       )

# Growth Stage look up table
GS_LUT <- tibble("GS" = c(0,1,2,3,4),
                 "GROWTH_STAGE" = c("Unknown","Juvenile","Adolescent","Mature","Old")
                 )

# EFG to TFI attributes look up table 
# read csv version of CGDL lookup table
TFI_LUT <- read.csv("./ReferenceTables/EFG_EVD_TFI.csv")[,c("EFG_NUM","MIN_LO_TFI",
                                                            "MIN_HI_TFI","MAX_TFI","EFG_NAME")]
names(TFI_LUT)[1] <- "EFG"
if (REGION_NO == 7)
  {clipPoly = adHocPolygon}else
  {clipPoly = "./ReferenceShapefiles/LF_DISTRICT.shp"}

inputR <- inputRasters(x = RasterRes)
outputFH <- file_path_sans_ext(basename(rawFH))


toc()#"initial setup")


## -------------------------------------------------------------------------
## Fire History processing -------------------------------------------------
print("Started Main FHanalysis")
tic("main FH analysis processing")

# Run Fire History analysis (function from EcoResFunctionsFMR)
FHanalysis <- FHProcess(rawFH = rawFH,
                        start.SEASON = start.SEASON,
                        end.SEASON = NULL,
                        OtherAndUnknown = OtherAndUnknown)                  #####-----use look up table if modified-----#####

# Add additional information to FHanalysis ???list/dataframe???             #####-----comment-----#####
FHanalysis$FireScenario = rawFH
FHanalysis$RasterRes = RasterRes
FHanalysis$ClipPolygonFile = clipPoly
FHanalysis$Region_No = REGION_NO
FHanalysis$PUBLIC_ONLY = PUBLIC_LAND_ONLY
FHanalysis$Start_Season = NULL
FHanalysis$name = paste0("FH_Analysis_", outputFH)

# Export Fire History Analysis as shapefile
st_write(FHanalysis$OutDF,
         file.path(ResultsDir, paste0(FHanalysis$name, ".shp"))
         )

toc()#"main FH analysis processing")

## Crop output rasters ------------------------------------------------------
tic("cropraster processing")
# crop rasters (function from EcoResFunctionsFMR)
cropRasters <- makeCropDetails(REG_NO = REGION_NO,
                               RasterRes = RasterRes,
                               PUBLIC_LAND_ONLY = PUBLIC_LAND_ONLY,
                               myPoly =clipPoly,
                               generalRasterDir = "./InputGeneralRasters"
                               )
# Add additional information to FHanalysis
FHanalysis$FH_IDr<-fasterize(sf=FHanalysis$OutDF,
                             raster =  cropRasters$Raster,
                             field = "ID",
                             fun="first"
                             )
save(FHanalysis,cropRasters,
     file = file.path(ResultsDir, paste0(FHanalysis$name,RasterRes,".rdata"))
     )

toc()#"cropraster processing")


## Combine all fire sequence data--------------------------------------------
# calculates and combines all fire sequence data from; FH_ID, EFG, and chosen polygons#####-----check------#####

#if reusing FH_processing from previous run: (details same as in settings file)#####-----delete------#####
#load(file.path(ResultsDir,paste0("FH_Analysis_",outputFH,RasterRes,".rdata")))#####-----delete------#####
#does "All combinations of fire sequences (FH_ID) and EFG plus chosen polygons #####-----delete------#####

tic("make all combs object")
# run function to calculate all combinations (function from calc_U_AllCombs)
myAllCombs <- calcU_All_Combs(FHAnalysis, cropRasters)

toc()#"make all combs object")


## Time between Fire Intervals (TFI) calculations ---------------------------
# TFI calculations (new calc_TFI_2)                                            #####-----can you version control this 'new', and remove _2/new reference------#####
tic("calc_TFI_2")
# run function to calculate TFI rasters (function from TFI_functionsFMR)
myTFI_2 <- calc_TFI_2(FHanalysis,
                      U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI,
                      Index_AllCombs = myAllCombs$Index_AllCombs,
                      TFI_LUT,
                      OutputRasters = makeTFIRasters)
fwrite(myTFI_2,
       file.path(ResultsDir, "TFI_Summary_2.csv"),
       na="",
       row.names = FALSE
       )

toc()#"calc_TFI_2")


## BBTFI calculations -----------------------------------------------------
# description                                                              #####-----comment-----#####
tic("BBTFI calculations complete")
# run function to calculate BBTFI (function from )
myBBTFI <- calcBBTFI_2(FHanalysis,
                       U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI,
                       Index_AllCombs = myAllCombs$Index_AllCombs,
                       TFI_LUT,
                       makeBBTFIrasters = makeBBTFIrasters)

# write BBTFI_LONG data
fwrite(myBBTFI$BBTFI_LONG,
       file.path(ResultsDir, "BBTFI_LONG.csv"),
       na = "",
       row.names = FALSE
       )

# write BBTFI_WIDE data
fwrite(myBBTFI$BBTFI_WIDE,
       file.path(ResultsDir, "BBTFI_WIDE.csv"),
       na = "",
       row.names = FALSE
       )

toc()#"BBTFI calculations complete")


## Gxxxx Sxxxx (GS) calculations ------------------------------------------
# description                                                              #####-----comment-----#####
tic("GS calculations")

#run function to calculate GS data (function from GS_Calcs)
GS_Summary <- makeGS_Sum(writeGSRasters,
                         ResultsDir,
                         U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI,
                         Index_AllCombs = myAllCombs$Index_AllCombs,
                         FHanalysis,
                         writeYears = NULL)

# write GS Summary VAT data
fwrite(GS_Summary$GS_Summary_wide,
       file.path(ResultsDir, "GS_Summary_VAT.csv"),
       na = "",
       row.names = FALSE
       )

# write GS Summary Long data
fwrite(GS_Summary$GS_Summary_Long,
       file.path(ResultsDir,"GS_Summary_Long.csv"),
       na="",
       row.names = FALSE
       )

toc()#"GS calculations complete")


## Garbage collection to free up memory -----------------------------------
gc()








#first sections of main FAME_FMR script to run with preexisting settings and FH_Analysis------------- 

# Clear workspace ---------------------------------------------------------


rm(list=ls(all=TRUE))
gc()
# Load libraries and custom functions -------------------------------------


options(stringsAsFactors = F)
library(Matrix.utils)
library(tools)
library (tidyr)
library (raster)
library(fasterize)
library(sf)
library(foreach)
library(doParallel)
library(dplyr)
library(tiff)
library(reshape2)
library(tabularaster)
library(tictoc)
library(filematrix)
library(data.table)
library(Rfast)
library(tictoc)
#loads functions used in TFI and RA calculations
source("EcoResFunctionsFMRv2.r")
source("TFI_functionsFMRv2.r")
source("calc_U_AllCombs.r")
source("GS_Calcs.R")
source("calcBBTFI_2.R")

#source("calc_TFI_2.r")
tic("whole process time")

#source("ButtonDisableHelpers.r")
#Set the maximum size of files for upload/ download 

#options(shiny.maxRequestSize=2*1024^3) 

# Load all settings and input files and tables ----------------------------
tic("whole process time")
source("./settings.r")


#MAKE RESULTS DIRECTORIES
#create a unique results directory for each run of scenarios 
#using starting time as a string for Results directory name

WD<-getwd()
#cleans up old resultsdirs which are empty
removeEmptyDirs(rootDir = "./")



for (d in c(ResultsDir)){dir.create(d)}
rm(d)
dir.create(file.path(ResultsDir,"RA_Rasters"))
dir.create(file.path(ResultsDir,"TFI_Rasters"))
dir.create(file.path(ResultsDir,"GS_Rasters"))

# Lookup table for choice of Fire region/ state or adhoc polygon for analysis----------

REG_LUT<-tibble(FIRE_REG = c(99, 1, 2, 3, 4, 5, 6, 7),
                FIRE_REGION_NAME = c("WHOLE OF STATE",
                                     "BARWON SOUTH WEST",
                                     "GIPPSLAND",
                                     "GRAMPIANS",
                                     "HUME",
                                     "LODDON MALLEE",
                                     "PORT PHILLIP",
                                     "USER DEFINED POLYGON")
)

#Lookup table for Fire FMZ-----------
FIREFMZ_LUT<-tibble(FIREFMZ = c(0, 1, 2, 3, 4, 5),
                    FIRE_FMZ_NAME = c("0 - Non FPA Land",
                                      "1 - Asset Protection Zone",
                                      "2 - Bushfire Moderation Zone",
                                      "3 - Landscape Management Zone",
                                      "4 - Planned Burn Exclusion Zone",
                                      "5 - Unknown. Contact Fire Management Officer"),
                    FIRE_FMZ_SHORT_NAME=c("Non FPA",
                                          "APZ",
                                          "BMZ",
                                          "LMZ",
                                          "PBEZ",
                                          "UNK")
)
#Lookup for delwp regeon names-----
DELWP_LUT<-tibble(DELWP= c(2,4,3,5,6,7),
                  DELWP_REGION=c("GIPPSLAND","HUME","PORT PHILLIP","BARWON SOUTH WEST","LODDON MALLEE","GRAMPIANS"))





#
#lookup table for TFI_STATUS values------------
TFI_STATUS_LUT<-structure(list(TFI_VAL = c(-99, 0, 1, 5, 6),
                               TFI_STATUS = c("NONE","WITHIN_TFI", "BELOW_MIN_TFI", "ABOVE_MAX_TFI", "ABOVE_MAX_BELOW_MIN_HI_TFI")),
                          row.names = c(NA,-5L), class = c("tbl_df", "tbl", "data.frame"))

write.csv(  TFI_STATUS_LUT,  file.path(ResultsDir,"TFI_Rasters","TFI_STATUS_LUT.csv"))

FIRETYPE_LUT<-tibble(  TYPE=c(1,2),  FIRETYPE=c("BURN","BUSHFIRE"))

#lookup for Growth Stage---------------
GS_LUT<-tibble("GS"=c(0,1,2,3,4),
               "GROWTH_STAGE"=c("Unknown","Juvenile","Adolescent","Mature","Old"))
#Lookup table from EFG to TFI attributes ( csv version of CGDL lookup table)
TFI_LUT<-read.csv("./ReferenceTables/EFG_EVD_TFI.csv")[,c("EFG_NUM","MIN_LO_TFI","MIN_HI_TFI","MAX_TFI","EFG_NAME")]
# ------------
names(TFI_LUT)[1]<-"EFG"
if (REGION_NO==7){clipPoly=adHocPolygon}else{clipPoly="./ReferenceShapefiles/LF_DISTRICT.shp"}

inputR<-inputRasters(x=RasterRes)
outputFH<-file_path_sans_ext(basename(rawFH))

load(file.path(ResultsDir,paste0("FH_Analysis_",outputFH,RasterRes,".rdata")))

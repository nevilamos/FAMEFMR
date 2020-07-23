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

tic("whole process time")
#source("ButtonDisableHelpers.r")
#Set the maximum size of files for upload/ download 

#options(shiny.maxRequestSize=2*1024^3) 

# Load all settings and input files and tables ----------------------------

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


# FH_processing starts here ------------------------------------------------

FHanalysis<-FHProcess(flattenedFH = rawFH,start.SEASON = start.SEASON,end.SEASON = NULL,OtherAndUnknown = OtherAndUnknown)
FHanalysis$FireScenario= rawFH
FHanalysis$RasterRes = RasterRes
FHanalysis$ClipPolygonFile = clipPoly
FHanalysis$Region_No = REGION_NO
FHanalysis$PUBLIC_ONLY = PUBLIC_LAND_ONLY
FHanalysis$Start_Season = NULL
FHanalysis$name<-paste0("FH_Analysis_",outputFH)
st_write(FHanalysis$OutDF,file.path(ResultsDir,paste0(FHanalysis$name,".shp")))


cropRasters<-makeCropDetails(REG_NO = REGION_NO,RasterRes = RasterRes,PUBLIC_LAND_ONLY = PUBLIC_LAND_ONLY,myPoly =clipPoly,generalRasterDir = "./InputGeneralRasters")
FHanalysis$FH_IDr<-fasterize(sf=FHanalysis$OutDF,raster =  cropRasters$Raster,field = "ID",fun="first")
save(FHanalysis,cropRasters,file=file.path(ResultsDir,paste0(FHanalysis$name,RasterRes,".rdata")))



# If reusing FH_processing from previous run ( details same as in settings file) -----------------------
#load(file.path(ResultsDir,paste0("FH_Analysis_",outputFH,RasterRes,".rdata")))
# TFI calcuations old calc_TFI ---------------------------------------------------------
tic("standard calcTFI")
myTFI<-calc_TFI(FHanalysis =FHanalysis,#the selected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                TFI_LUT_DF = TFI_LUT,
                cropRasters = cropRasters,
                OutputRasters = writeSpRasters)
print("Finished my TFI")
write.csv(myTFI,file.path(ResultsDir,"TFI_Summary.csv"))
toc()

# TFI calcuations new calc_TFI_2 ---------------------------------------------------------
tic("calc_TFI_2")
myTFI_2<-calc_TFI_2(FHanalysis =FHanalysis,#the selected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                TFI_LUT,
                cropRasters = cropRasters,
                OutputRasters = "No")
write.csv(myTFI_2,file.path(ResultsDir,"TFI_Summary_2.csv"))
toc()

# myTFI_lorequ<-calc_TFI_lorequ(FHanalysis =FHanalysis,#the selected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
#                 TFI_LUT_DF = TFI_LUT,
#                 cropRasters = cropRasters,
#                 OutputRasters = writeSpRasters)
# print("Finished my TFI")
# write.csv(myTFI_lorequ,file.path(ResultsDir,"TFI_lorequ_Summary.csv"))

# BBTFI calculations ------------------------------------------------------
tic("BBTFI calculations complete")
myBBTFI<-calc_BBTFI(FHanalysis,
                    cropRasters,#,the selected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                    TFI_LUT_DF = TFI_LUT)
print("finished BBTFI calcs")
save(myBBTFI,file=file.path(ResultsDir,paste(file_path_sans_ext(FHanalysis$name),"BBTFI_TFI.rdata")))
# 
# write.csv(myTFI$UNDER_TFI__BY_EFG_WIDE,file.path(ResultsDir,"UnderTFIbyEFGandSEASONwide.csv"),row.names=F)
# 
# toc()
# 
save(myBBTFI,file=file.path(ResultsDir,"bbtfi_tfi.rdata"))
write.csv(myBBTFI$BBTFI_BY_TYPE,file.path(ResultsDir,"BBTFI_BY_TYPE.csv"))
BBTFI_time<-toc()
print(BBTFI_time)
# GS calculations ---------------------------------------------------------




tic("GS calculations")
GS_LU<-makeGS_LU()


tsf_ysf_mat<-makeYSF_LFT_matrix(FHanalysis = FHanalysis,
                                myCropDetails=cropRasters,
                                FH_ID.tif=FHanalysis$FH_IDr)
gc()
GS_Summary<-makeGS_Sum(TimeSpan = FHanalysis$TimeSpan,
                       writeGSRasters = "Yes",
                       myLU = GS_LU,
                       myResultsDir = ResultsDir,
                       myCropDetails = cropRasters,
                       myFHResults = FHanalysis,
                       myYSF_LFT = tsf_ysf_mat,
                       writeYears=NULL)
write.csv(GS_Summary,file.path(ResultsDir,"GS_Summary.csv"))
gc()
toc()
TimeTaken<-toc()
#save.image(file.path(ResultsDir,paste0(myBasename,"_finalImage.rdata")))






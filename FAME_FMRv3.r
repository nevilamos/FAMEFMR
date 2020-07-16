
# Clear workspace ---------------------------------------------------------


rm(list=ls(all=TRUE))

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
#library(filematrix)
library(data.table)
library(Rfast)
#loads functions used in TFI and RA calculations
source("EcoResFunctionsFMRv2.r")
source("TFI_functionsFMRv2.r")
source("calc_TFI_75.r")
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

# Lookup table for choice of region/ state or adhoc plygon for analysis

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

#Lookup table for Fire FMZ
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

#lookup table for TFI_STATUS values
TFI_STATUS_LUT<-tibble(
  TFI_VAL= c(-99,0,1,2),
  TFI_STATUS=c("NA","WITHIN_TFI","BELOW_MIN_TFI","ABOVE_MAX_TFI")
)

write.csv(
  TFI_STATUS_LUT,
  file.path(ResultsDir,"TFI_Rasters","TFI_STATUS_LUT.csv"))

FIRETYPE_LUT<-tibble(
  TYPE=c(1,2),
  FIRETYPE=c("BURN","BUSHFIRE")
)
#lookup for Growth Stage
GS_LUT<-tibble("GS"=c(0,1,2,3,4),
               "GROWTH_STAGE"=c("Unknown","Juvenile","Adolescent","Mature","Old"))
#Lookup table from EFG to TFI attributes ( csv version of CGDL lookup table)
TFI_LUT<-read.csv("./ReferenceTables/EFG_EVD_TFI.csv")[,c("EFG_NUM","MIN_LO_TFI","MIN_HI_TFI","MAX_TFI","EFG_NAME")]

names(TFI_LUT)[1]<-"EFG"
if (REGION_NO==7){clipPoly=adHocPolygon}else{clipPoly="./ReferenceShapefiles/LF_DISTRICT.shp"}


#e = local({load("/mnt/TbSSD/FAMEStatewideReporting2020/ResultsFFRAU_FH2020_FAME_vg94_STATEWIDE_NO_1755_75m/FH_Analysis_FFRAU_FH2020_FAME_vg94_STATEWIDE_SEASON_AS_NUMERIC_NO_175575.rdata"); environment()})
#tools:::makeLazyLoadDB(e, "New")






# FH_procssing starts here ------------------------------------------------


#Starts Processing here
inputR<-inputRasters(x=RasterRes)
outputFH<-file_path_sans_ext(basename(rawFH))

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
#load(file.path(ResultsDir,paste0(FHanalysis$name,RasterRes,".rdata")))


# TFI calcuations ---------------------------------------------------------




myBasename<-file_path_sans_ext(basename(FHanalysis$name))


myTFI<-calc_TFI(FHanalysis =FHanalysis,#the selected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                TFI_LUT_DF = TFI_LUT,
                cropRasters = cropRasters,
                OutputRasters = writeSpRasters)
print("Finished my TFI")
write.csv(myTFI,file.path(ResultsDir,"TFI_Summary.csv"))

myTFIv75<-calc_TFI_75(FHanalysis =FHanalysis,#the selected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                TFI_LUT_DF = TFI_LUT,
                cropRasters = cropRasters,
                OutputRasters = writeSpRasters)
print("Finished my TFI")
write.csv(myTFIv75,file.path(ResultsDir,"TFI_Summary_v75.csv"))

# myTFI_lorequ<-calc_TFI_lorequ(FHanalysis =FHanalysis,#the selected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
#                 TFI_LUT_DF = TFI_LUT,
#                 cropRasters = cropRasters,
#                 OutputRasters = writeSpRasters)
# print("Finished my TFI")
# write.csv(myTFI_lorequ,file.path(ResultsDir,"TFI_lorequ_Summary.csv"))

# BBTFI calculations ------------------------------------------------------

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
#save(myBBTFI,file="bbtfi_tfi.rdata")


# GS calculations ---------------------------------------------------------


write.csv(myBBTFI$BBTFI_BY_TYPE,file.path(ResultsDir,"BBTFI_BY_TYPE.csv"))
GS_LU<-makeGS_LU()


tsf_ysf_mat<-makeYSF_LFT_matrix(FHanalysis = FHanalysis,
                                myCropDetails=cropRasters,
                                FH_ID.tif=FHanalysis$FH_IDr)
gc()
GS_Summary<-makeGS_Sum(TimeSpan = FHanalysis$TimeSpan,
                       #myHDMSpp_NO = HDMSpp_NO,
                       writeGSRasters = "Yes",
                       myLU = GS_LU,
                       myResultsDir = ResultsDir,
                       myCropDetails = cropRasters,
                       myFHResults = FHanalysis,
                       myYSF_LFT = tsf_ysf_mat,
                       writeYears=NULL)
write.csv(GS_Summary,file.path(ResultsDir,"GS_Summary.csv"))
gc()

save.image(file.path(ResultsDir,paste0(myBasename,"_finalImage.rdata")))
#########Species based stuff
# LU_List<-makeLU_List(myHDMSpp_NO = HDMSpp_NO,
#myAbundDataLong = AbundDataLong)
# if(is.null(customSpList)){
#   TaxonList <-read.csv("./ReferenceTables/DraftTaxonListStatewidev2.csv")
# }else{
#   TaxonList<-read.csv(customSpList)
# }
# HDMSpp_NO<<-TaxonList$TAXON_ID[TaxonList$Include=="Yes"]
# 
# print("getting HDMvals")
# if(FHanalysis$RasterRes==225){
#   load(paste0("./HDMS/HDMVals",RasterRes,".rdata"))
#   HDMVals<-HDMVals[cropRasters$IDX,as.character(HDMSpp_NO)]
# }else{
#   print("doing 75m version makeHDMValsfromRasters")
#   HDMVals<<-makeHDMValsfromRasters(myHDMSpp_NO = HDMSpp_NO,
#                                    myCropDetails = cropRasters)
# }
# 
# 
# if(is.null(customResponseFile)){
#   mySpGSResponses="./ReferenceTables/OrdinalExpertLong.csv"
# }else{
#   mySpGSResponses=customResponseFile
# }
# AbundDataByGS = read.csv(mySpGSResponses)[,c("EFG_NO", "GS4_NO",  "FireType" , "Abund", "VBA_CODE")]  #Select the file giving the fauna relative abundacne inputs you wish to use
# AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="High"]<-2
# AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="Low"]<-1
# AbundDataByGS<-AbundDataByGS[!is.na(AbundDataByGS$Abund),c("EFG_NO", "GS4_NO",  "FireTypeNo" , "Abund", "VBA_CODE")]
# 
# EFG_TSF_4GS <- read.csv("./ReferenceTables/EFG_TSF_4GScorrectedAllEFGto400yrs.csv")[,c('EFG_NO','GS4_NO',"YSF")]
# AbundDataLong = merge(AbundDataByGS, EFG_TSF_4GS,   by=c('EFG_NO','GS4_NO'))
# AbundDataLong<-AbundDataLong[order(AbundDataLong$VBA_CODE),]

# SpYearSumm<-makeSppYearSum(TimeSpan = FHanalysis$TimeSpan,
#                            myHDMSpp_NO = HDMSpp_NO,
#                            writeSpRasters = writeSpRasters,
#                            myLU_List = LU_List,
#                            YSF_TSF_Dir = YSF_TSF_Dir,
#                            ResultsDir = ResultsDir,
#                            EFG = cropRasters$EFG,
#                            myCropDetails = cropRasters,
#                            HDMVals = HDMVals,
#                            myFHResults = FHanalysis,
#                            myYSF_LFT = tsf_ysf_mat,
#                            TaxonList = TaxonList,
#                            writeYears = yearsForRasters)
# write.csv(SpYearSumm,file.path(ResultsDir,"SpYearSumm.csv"),row.names=F)
# myBaseline<<-ifelse(endBaseline>startBaseline,startBaseline:endBaseline,input$startBaseline)
# calcDeltaAbund(SpYearSumm,
#                TimeSpan=FHanalysis$TimeSpan,
#                myBaseline,
#                ResultsDir,
#                HDMSpp_NO,
#                TaxonList)
# makeSummaryGraphs(SpYearSumm,
#                   ResultsDir,
#                   HDMSpp_NO,
#                   outputFH<-outputFH)
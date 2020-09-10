#this file is the script that runs all the functions concerned with calcuation species post fire abundances


# Choice of complete or custom species list determining which species to include in analysis
# if a custom species list is used this can only contain species for which there is a raster in the HDMS files and species for which there 
# this file must follow the formatting of the deafualt species list 
#is a non- zero abundance in the extent of the analysis
#If rasters are added to or updated these files then the process to update the HDMValue filematrices must be run run using the script "makeHDM_Value_Matrices.r"
tic()
if(is.null(customSpList)){
  TaxonList <-read.csv("./ReferenceTables/DraftTaxonListStatewidev2.csv")#default species list
}else{
  TaxonList<-read.csv(customSpList)#custom species list
}

#filters the species to be considered in the analysis - these are flagged with the value "Yes" in the incude column in the input species list.
HDMSpp_NO<-TaxonList$TAXON_ID[TaxonList$Include=="Yes"]
#Gets a subset of those for which rasters are to be written when WriteSpeciesRasters =="Yes" If this is NULL then all are written
writeSp<-TaxonList$TAXON_ID[TaxonList$WriteSpeciesRaster=="Yes"]
writeSp<-writeSp[writeSp%in%HDMSpp_NO]
tic("time to read in species values")
print("getting HDMvals")
if(FHanalysis$RasterRes==225){
  load(paste0("./HDMS/HDMVals",RasterRes,".rdata"))
  HDMVals<-HDMVals[cropRasters$IDX,as.character(HDMSpp_NO)]
}else{
  print("doing 75m version makeHDMValsfromRasters")
  HDMVals<<-makeHDMValsfromRasters(myHDMSpp_NO = HDMSpp_NO,
                                   myCropDetails = cropRasters)
}
toc()
# choice of default  or custom species abundance values input files
if(is.null(customResponseFile)){
  mySpGSResponses="./ReferenceTables/OrdinalExpertLong.csv" #default valeus ( from FFO database)
}else{
  mySpGSResponses=customResponseFile # alternative user define custom abundance valuesin matching format
}
tic("making AbundDataLong")
#the below lines inflate the input file containing species abundance by growth stage to a data frame (AbundDataLong) containing specieabundance (still stepped by growth stage)
###this sections would be skipped if a file for species response by time since fire ( rather than growth stage is provided ############
AbundDataByGS = read.csv(mySpGSResponses)[,c("EFG_NO", "GS4_NO",  "FireType" , "Abund", "VBA_CODE")]  #Select the file giving the fauna relative abundance inputs you wish to use
AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="High"]<-2
AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="Low"]<-1
AbundDataByGS<-AbundDataByGS[!is.na(AbundDataByGS$Abund),c("EFG_NO", "GS4_NO",  "FireTypeNo" , "Abund", "VBA_CODE")]
EFG_TSF_4GS <- read.csv("./ReferenceTables/EFG_TSF_4GScorrectedAllEFGto400yrsV2.csv")[,c('EFG_NO','GS4_NO',"YSF")]
AbundDataLong = merge(AbundDataByGS, EFG_TSF_4GS,   by=c('EFG_NO','GS4_NO'))

###if species response by time since fire file is provided uncomment following line and comment out section above back to last comment.
#This file must go u to 401 years since fire for all species and follow the formatting of AbundDatLong object exactly.
#the file will need to be defined in settings.

#AbundDataLong<-read.csv(customAbundanceByTSF_EFG_FT_file)

AbundDataLong<-AbundDataLong[order(AbundDataLong$VBA_CODE),]
toc()
tic("making Spp abund LU List")
LU_List<-makeLU_List(myHDMSpp_NO = HDMSpp_NO,
                      myAbundDataLong = AbundDataLong)
toc()

tic("main species abundace calulations in makeSppYearSumm2")
SpYearSummWide<-makeSppYearSum2(FHanalysis,
                                      myHDMSpp_NO = HDMSpp_NO,
                                      writeSpRasters = writeSpRasters,
                                      myLU_List = LU_List,
                                      ResultsDir = ResultsDir,
                                      HDMVals = HDMVals,
                                      TaxonList = TaxonList,
                                      writeYears=NULL,
                                      writeSp =writeSp
)
toc()
myBaseline<-ifelse(endBaseline>startBaseline,startBaseline:endBaseline,startBaseline)
tic("calculate changes in relative abundance")
myDeltaAbund<-calcDeltaAbund(SpYearSummSpreadbyYear =SpYearSummWide,
                          TimeSpan =FHanalysis$TimeSpan,
                          myBaseline,
                          ResultsDir,
                          HDMSpp_NO,
                          TaxonList)
toc()
toc()


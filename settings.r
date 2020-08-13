#R Settings file for FAMEv1.0.1 batch process. Only the section after the "=" in line 
# 8 + should be edited.  the file should be saved as settings.r  in the input directory, all other input file names required should be in qutes, and also in the input directory
# At the moment this does not run the aspatial GSO tool.  it runs all the spatial TFI and GSO calcuations. if customSppList and customAbundanceLU =NULL then the default files are used. 
# if REGION_NO =7 then an adHocPolygon shapefile must also be provided ( again in the inputs directory.)


#the integer value of the region number 1-6 for FFR regions,7 for user suppied adHoc polygon 99 for Statewide
REGION_NO 			= 99 
#path to the ad Hoc polygon if REGION_NO 			== 7
#adHocPolygon 		= "./rawFH/BBTFI_checkArea.shp"
#path to your rawFH file ( output of the ARCGIS preprocessing tool)
  rawFH 				= "./rawFH/FFRAU_FH2020_FAME_vg94Draft_JFMP_for_Eco_Assessment_Burnabclip2022_SEAS_0_Removed.shp"
#path to custom species list if used
customSpList 		= NULL
#path to custom response list if used
customResponseFile 	= NULL
#raster resolution 75 or 225
RasterRes			= 225
#sets the csv table containing 4GS info by years since fire
myEFG_TSF_4GS= read.csv("./ReferenceTables/EFG_TSF_4GScorrectedAllEFGto400yrsV2.csv")[,c('EFG_NO','GS4_NO',"YSF")]
#path to HDMVals225 file
HDMVals225 			=	"./HDMS/HDMVals225.rdata"
#whether to output species rasters
writeSpRasters		= "No" # "No"
#whether to output species rasters
writeGSRasters		= "No" # "No"
#if writeSpRasters		= "Yes" vector of years for which rasters are to be written
yearsForRasters = NULL
#whether to write TFI rasters
makeTFIRasters = "N0"#Yes"

#whether to output bbtfi rasters
makeBBTFIrasters=TRUE#FALSE
#"first season for which output is wanted ( four digit year as integer)
#if NUll then second season in in history  is used cannot use first season because it has no interval, this may still fail if there is no overlap,
start.SEASON		= 1980 
#value to use for cases where fire type is "OTHER" or "UNKNOWN",1 ="BURN",2="BUSHFIRE",NA = Fire excluded from analysis default is 2 ("BUSHFIRE")
OtherAndUnknown 	= 2# (2,1,NA)
#whether analysis should be undetaken only on public land
PUBLIC_LAND_ONLY	=	FALSE
#start and end baseline years if single year then the two values should be equal
startBaseline =1990
endBaseline =1990

Ncores=4
print(paste("Using",Ncores,"cores"))

  ResultsDir<-file.path("./Results",paste(tools::file_path_sans_ext(basename(rawFH)),RasterRes,sep="_"))
  dir.create(ResultsDir)
file.copy("settings.r",file.path(ResultsDir,"settings.r"))
#change beleow to TRUE for abrreviated test ( on first two years only)
Test=F

source("global.r")
unlink(resultsDir,recursive = TRUE)
#source("batchsettings_JFMP.r")
source("batchsettings_raster_first_fh.r")


#(rv name in kept for consistency with rv <- reactiveValues in shiny app)
# LOAD settings ----
rv<-list()
#  rawFH file to be run ----
rv$rawFHPath <- rawFHPaths


rv$rawFHName <- basename(rv$rawFHPath)

#AdHoc shapefile file to be run ----
rv$AdHocPath <- AdHocPath
if(!is.null(rv$AdHocPath)){
  rv$AdHocName <- basename(rv$AdHocPath)
}
# whether or not to make make TFI rasters----
rv$makeTFIrasters <- makeTFIrasters

# whether or not to make  BBTFI rasters----
rv$makeBBTFIrasters <- makeBBTFIrasters

# whether or not to make  growth stage  rasters----
rv$makeGSrasters <- makeGSrasters

#whether or not to run species calculations

#customSpList be run ----

rv$customSpList <- customSpList
rv$customSpListName <- basename(rv$customSpList)

# Observer of choice for PU polys sets values to null if FALSE  ----
rv$usePUpolys = usePUpolys
if (!rv$usePUpolys == TRUE) {
  rv$puPath <- NULL
  rv$puName <- NULL
} else {
  # PU shapefile file to be run ----
  rv$puPath <- puPath
  rv$puName <- basename(rv$puPath)
  if (runJFMP == TRUE){
    # zone wt File be run ----
    rv$zoneWtFile <- zoneWtFile
    rv$zoneWtFileName <- basename(rv$zoneWtFile)
    rv$zoneWt <- read_csv(rv$zoneWtFile)
    # jfmp metric wt File be run ----
    rv$jfmpMetricWtFilePath <- jfmpMetricWtFile
    rv$jfmpMetricWtFileName <- basename(rv$jfmpMetricWtFilePath)
    rv$jfmpMetricWt <- read_csv(rv$jfmpMetricWtFilePath)
    # JFMP Area Target file  ----

    rv$targetHaFilepath <- targetHaFilepath
    rv$targetHaFileName <- basename(rv$targetHaFilepath)
    rv$targetHa <- read_csv(rv$targetHaFilepath)
  }

}

rv$spResponseChoice <- spResponseChoice

if (rv$spResponseChoice == TRUE){
  rv$customResponseFile <- customResponseFile
  rv$customResponseName <- basename(rv$customResponseFile)
}
rv$spListChoice <- spListChoice
rv$spResponseChoice <- spResponseChoice
rv$abundByGS <- abundByGS
rv$makeRArasters <- makeRArasters

# NUMERIC SETTINGS  ----
rv$RasterRes <- RasterRes
rv$startTimespan <- startTimespan
rv$endTimespan <- endTimespan

rv$startBaseline <- startBaseline

rv$endBaseline <- endBaseline

rv$public <- public

rv$otherAndUnknown <- otherAndUnknown

# SEASONS for write rasters (yearsForRasters)  ----

rv$yearsForRasters <- yearsForRasters

rv$REGION_NO <- REGION_NO


rv$outputFH <- file_path_sans_ext(basename(rv$rawFHPath))
if (rv$usePUpolys == TRUE) {
  rv$outputFH <- paste(rv$outputFH, file_path_sans_ext(basename(rv$puName)), sep = "_")
}

#output directories creation

rv$resultsDir<-file.path(resultsDir,rv$outputFH)
for (i in c("RA_Rasters","TFI_Rasters","BBTFI_Rasters")){
  dir.create(file.path(rv$resultsDir,i),recursive = TRUE)
}
myREG_NO <- as.integer(rv$REGION_NO)
RasterRes <- as.integer(rv$RasterRes)
print(paste("RasterRes =", RasterRes))
HDM_RASTER_PATH <-
  fame_data_path("HDMS", paste0(rv$RasterRes, "m/BinaryThresholded"))

if (rv$REGION_NO == 7) {
  rv$clipShape <- rv$AdHocPath
} else {
  rv$clipShape <- fame_data_path("ReferenceShapefiles","LF_REGION.shp")
}


rv$end.SEASON <- rv$endTimespan


#alteration to settings if doing JFMP1
if (runJFMP  == TRUE){
  rv$JFMPSeason0= JFMPSeason0
  rv$end.SEASON <- rv$JFMPSeason0 + 4
  rv$makeRArasters = FALSE
  rv$makeGSrasters = FALSE
  rv$makeTFIrasters = FALSE
  doSpeciesCalculations = FALSE
  rv$spListChoice =TRUE
  rv$TaxonListPath = jfmpSpeciesList
  rv$customSpList <- jfmpSpeciesList
  rv$customSpListName <- basename(rv$customSpList)
}


#Starts analysis here ----

#make cropRAsters list of extent and values for rasters cropped to extent of area of interest. -----
rv$cropRasters <- cropToOutput(
  REG_NO = rv$REGION_NO,
  myRasterRes = rv$RasterRes,
  PUBLIC_LAND_ONLY = rv$public,
  myPoly = rv$clipShape,
  generalRasterDir = fame_data_path("InputGeneralRasters")
)

rv$cropRasters$HDM_RASTER_PATH <- HDM_RASTER_PATH
template_r<-eval(rv$cropRasters$rasterDef)
values(template_r)<-rv$cropRasters$RGN>0
# main Fire history process ------
rv$FHAnalysis <- fhProcess_raster_first(
  r_template = template_r,
  input_fire_history = rv$rawFHPath,
  start.SEASON = rv$startTimespan,
  end.SEASON = rv$end.SEASON,
  OtherAndUnknown = rv$otherAndUnknown
)
# Append into FH analysis object
rv$FHAnalysis$FireScenario <- rv$rawFHName
rv$FHAnalysis$RasterRes <- rv$RasterRes
rv$FHAnalysis$ClipPolygonFile <- rv$clipShape
rv$FHAnalysis$Region_No <- myREG_NO
rv$FHAnalysis$PUBLIC_ONLY <- rv$public
rv$FHAnalysis$Start_Season <- NULL

rv$FHAnalysis$name <- paste0("FH_Analysis_", rv$outputFH)

st_write(rv$FHAnalysis$OutDF,
         file.path(rv$resultsDir, paste0(rv$FHAnalysis$name, ".",tools::file_ext(rv$FHAnalysis$FireScenario))),
         append = FALSE
)


rv$cropRasters$FH_ID<-rv$FHAnalysis$FH_ID
#rv$cropRasters$FH_ID<-values(rv$FHAnalysis$FH_IDr)
#mask  region values ( which is masked by plm where public land only)
rv$cropRasters$FH_ID[is.na(rv$cropRasters$RGN)]<-NA


# check if pupoly is to be used
if (rv$usePUpolys == TRUE){
  myPuPoly <- rv$puPath
  if (runJFMP==T) {
    # update the FHAnalysis$OutDF  with noburn columns
    rv$FHAnalysis$OutDF <-
      rv$FHAnalysis$OutDF %>% bind_cols(make_JFMPNoBurnTab(
        myFHAnalysis = rv$FHAnalysis,
        JFMPSeason0 = rv$JFMPSeason0
      ))
    rv$FHAnalysis$YSFNames <- c(rv$FHAnalysis$YSFNames, "YSFNoBurn")
    rv$FHAnalysis$LBYNames <- c(rv$FHAnalysis$LBYNames, "LBYNoBurn")
    rv$FHAnalysis$LFTNames <- c(rv$FHAnalysis$LFTNames, "LFTNoBurn")

    print("appended JFMPNoBurnCols")
  } }else {
    myPuPoly <- NULL

  }




rv$allCombs <- calcU_All_Combs(
  myFHAnalysis = rv$FHAnalysis,
  myCropRasters = rv$cropRasters,
  myRasterRes = RasterRes,
  puPoly = rv$puPath
)

print("made allcombs")

r<-eval(rv$cropRasters$rasterDef)


print("finished FH analysis")
#save completed FHanalysis rv with
analysisPath<-file.path(rv$resultsDir ,paste0(rv$outputFH,"_FHanalysis",rv$RasterRes,".qs"))
saveSpatRasterList(rv, analysisPath)


#block for TFI and GS calcuations in try()----

# TFI  related calculations----

print("running TFI calc")

rv$TFI <- calc_TFI_2(
  myFHAnalysis = rv$FHAnalysis,
  myCropRasters = rv$cropRasters,
  myAllCombs = rv$allCombs,
  myTFI_LUT = TFI_LUT,
  OutputRasters = T,
  myResultsDir = rv$resultsDir
)



# need to change the sort order for the factor to get correct stacking
# order ( no alphabetical) on chart
rv$TFI$TFI_STATUS <-
  factor(
    rv$TFI$TFI_STATUS,
    levels = c(
      "BELOW_MIN_TFI",
      "WITHIN_TFI",
      "ABOVE_MAX_TFI",
      "NONE"
    )
  )



# write results out to csv files
readr::write_csv(rv$TFI,
                 file = file.path(rv$resultsDir, "TFI_LONG.csv")
)
readr::write_csv(
  rv$TFI %>%
    group_by(EFG_NAME, SEASON, TFI_STATUS) %>%
    summarise(AreaHa = sum(Hectares)) %>%
    pivot_wider(
      names_from = SEASON,
      values_from = AreaHa,
      values_fill = 0
    ),
  file = file.path(rv$resultsDir, "TFI_EFG_SUMMARY.csv")
)

print("Finished TFI calculations")

# BBTFI  related calculations----


print("calculating BBTFI")
#source("D:/FAMEFMR/R/calcBBTFI_2.R")
rv$BBTFI <- calcBBTFI_2(
  myFHAnalysis = rv$FHAnalysis,
  myAllCombs = rv$allCombs,
  myCropRasters = rv$cropRasters,
  makeBBTFIrasters = rv$makeBBTFIrasters,
  myResultsDir = rv$resultsDir
)
print("finished BBTFI calcs")



write.csv(rv$BBTFI$BBTFI_LONG,
          file = file.path(
            rv$resultsDir,
            "BBTFI_LONG.csv"
          ),
          row.names = FALSE
)
write.csv(
  rv$BBTFI$BBTFI_LONG %>%
    group_by(EFG_NAME, TBTFI) %>%
    summarise(AreaHa = sum(Hectares)) %>%
    pivot_wider(names_from = TBTFI, values_from = AreaHa),
  file = file.path(
    rv$resultsDir,
    "TimesBBTFI_SUMMARY.csv"
  ),
  row.names = FALSE
)

write.csv(rv$BBTFI$BBTFI_WIDE,
          file = file.path(
            rv$resultsDir,
            "BBTFI_WIDE.csv"
          ),
          row.names = FALSE
)

# Run GS calculations----


print("GS Calculations")
GS_Summary<-calc_growth_stages(myFHAnalysis = rv$FHAnalysis,
                               myCropRasters = rv$cropRasters,
                               myAllCombs = rv$allCombs,
                               OutputRasters = rv$makeGSrasters,
                               myResultsDir = rv$resultsDir)


rv$GS_Summary <- GS_Summary

readr::write_csv(rv$GS_Summary$GS_Summary_Long,
                 file = file.path(
                   rv$resultsDir,
                   "GS_LONG.csv"
                 )
)

readr::write_csv(rv$GS_Summary$GS_Summary_wide,
                 file = file.path(
                   rv$resultsDir,
                   "GS_WIDE.csv"
                 )
)
print("finished GS calcs")


# save output so far----
pathToDelete<-analysisPath
analysisPath<-file.path(rv$resultsDir ,paste0(rv$outputFH,"_TFI_GS",".qs"))
saveSpatRasterList(rv, analysisPath)
unlink(pathToDelete)



# species calculations ------
if(doSpeciesCalculations){
  rv<-run_relative_abundance_analysis()
  }
if(runJFMP==TRUE){
  rv<-run_relative_abundance_analysis()
  validate(
    need(rv$usePUpolys == TRUE, message = "The FHAnalysis does not contain planning/burn units for JFMP calculations")
  )
  validate(need(length(rv$jfmpMetricWtFilePath) > 0, message = "You need to select a JFMP metric weight file"))
  validate(need(length(rv$zoneWtFile) > 0, message = "You need to select a JFMP zone Weight file"))
  print("doing JFMP1")
  # read in lookup tables for weighting of JFMP
  rv$zoneWt <- read_csv(rv$zoneWtFile)
  rv$jfmpMetricWt <- read_csv(rv$jfmpMetricWtFilePath)

  puDF <- jfmp1(
    myPUPath = rv$puPath,
    grpSpYearSumm = rv$SpYearSumm$grpSpYearSumm,
    myAllCombs = rv$allCombs,
    myTaxonList = rv$TaxonList,
    myBBTFI = rv$BBTFI,
    myJFMPSeason0 = rv$JFMPSeason0,
    zoneWt = rv$zoneWt,
    jfmpMetricWt = rv$jfmpMetricWt
  )
  rv$puDF <- puDF
  readr::write_csv(
    rv$puDF,
    file.path(
      rv$resultsDir,
      paste0(
        "PU_Rankings_",
        tools::file_path_sans_ext(rv$puName),
        ".csv"
      )
    )
  )
  gc()

  print(" jfmp1 PU_Rankings output saved to disk")

  print("running autoJFMP")
  rv$targetHa <- read_csv(rv$targetHaFilepath)
  rv$autoJFMP <- autoJFMP(
    myJFMP1 = rv$puDF,
    myTargetHa = rv$targetHa
  )

  readr::write_csv(
    rv$autoJFMP,
    file.path(
      rv$resultsDir,
      paste0(
        "AutoJFMP_",
        tools::file_path_sans_ext(rv$puName),
        ".csv"
      )
    )
  )

  autoJFMPsummary <- jfmpSummary(rv$autoJFMP)
  rv$autoJFMPsummary <- autoJFMPsummary

  readr::write_csv(
    rv$autoJFMPsummary,
    file.path(
      rv$resultsDir,
      paste0(
        "AutoJFMP_summary_",
        tools::file_path_sans_ext(rv$puName),
        ".csv"
      )
    )
  )

  print("finished auto JFMP")



}






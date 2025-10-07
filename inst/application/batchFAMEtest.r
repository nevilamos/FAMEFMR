source("global.r")
unlink(resultsDir,recursive = TRUE)
source("batchsettings.r")


#(rv name in kept for consistency with rv <- reactiveValues in shiny app)
# LOAD settings ----
rv<-list()

  #  rawFH file to be run ----
  rv$rawFHPath <- firstFH
  rv$firstFHLayer <-firstFHLayer
  rv$secondFH=secondFH
  rv$secondFHLayer=secondFHLayer
  rv$rawFHName <- basename(rv$rawFHPath)
  rv$baseFire <-baseFire

  
  #AdHoc shapefile file to be run ----
  rv$AdHocPath <- AdHocPath
  if(!is.null(rv$AdHocPath)){
    rv$AdHocName <- basename(rv$AdHocPath)
  }
  # whether or not to make make TFI rasters----
  rv$makeTFIrasters <- makeTFIrasters
  
  # whether or not to make  BBTFI rasters----
  rv$makeBBTFIrasters <- makeBBTFIrasters
  
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
    rv$jfmpMetricWtFile <- jfmpMetricWtFile
    rv$jfmpMetricWtFileName <- basename(rv$jfmpMetricWtFile)
    rv$jfmpMetricWt <- read_csv(rv$jfmpMetricWtFile)
    # JFMP Area Target file  ----
    
    rv$targetHaFilepath <- targetHaFilepath
    rv$targetHaFileName <- basename(rv$targetHaFilepath)
    rv$targetHa <- read_csv(rv$targetHaFilepath)
    }

  }
  
  # Observer of custom relative abundance table choice----
  
  rv$spResponseChoice <- spResponseChoice
  
  if (rv$spResponseChoice == TRUE){
    # customResponseFile be run ----
    
    rv$customResponseFile <- customResponseFile
    rv$customResponseName <- basename(rv$customResponseFile) 
  }
  
  
  
  
  # Observer of custom species list choice  ----
  
  rv$spListChoice <- spListChoice
  
  
  
  # Observer of custom relative abundance table choice----
  
  rv$spResponseChoice <- spResponseChoice
  
  
  # Observer of Relative abundance table by growth stage choice ----
  
  rv$abundByGS <- abundByGS
  
  # Observer of make rasters choice ----
  rv$makeRArasters <- makeRArasters
  
  
  
  
  # OBSERVERS of NUMERIC SETTINGS  ----
  # Observer for RasterRes----
  
  rv$RasterRes <- RasterRes
  
  
  # Observer for First season for analysis output (startTimespan)----
  
  rv$startTimespan <- startTimespan
  rv$endTimespan <- endTimespan
  # Observer for start baseline----
  rv$startBaseline <- startBaseline
  
  # Observer for end baseline----
  
  rv$endBaseline <- endBaseline

  # set public land----
  
  rv$public <- public

  rv$otherUnknown <- otherUnknown

  # SEASONS for write rasters (yearsForRasters)  ----
  # this is not working at the moment to reload the seasons selected in previous session
  
  rv$yearsForRasters <- yearsForRasters
  

  #choose a region  ----
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
    paste0("./HDMS/", rv$RasterRes, "m/BinaryThresholded")
  
  if (rv$REGION_NO == 7) {
    rv$clipShape <- rv$AdHocPath
  } else {
    rv$clipShape <- "./ReferenceShapefiles/LF_DISTRICT.shp"
  }
  
  
  if (runJFMP == TRUE) {
    rv$end.SEASON <- rv$JFMPSeason0 + 4
  } else {
    rv$end.SEASON <- rv$endTimespan
  }
  
  rv$max_interval<-max_interval
#Starts analysis here ---- 
  
#make cropRAsters list of extent and values for rasters cropped to extent of area of interest. ----- 
  rv$cropRasters <- cropToOutput(
    REG_NO = myREG_NO,
    myRasterRes = RasterRes,
    PUBLIC_LAND_ONLY = rv$public,
    myPoly = rv$clipShape,
    generalRasterDir = "./InputGeneralRasters"
  )
  
  rv$cropRasters$HDM_RASTER_PATH <- HDM_RASTER_PATH

  rv$FHAnalysis<- fhProcess2(fhProcess1(inFH = rv$rawFHPath,
                                        inFHLayer = rv$inFHLayer,
               OtherAndUnknown = 2,
               validFIRETYPE = c("BURN", "BUSHFIRE", "UNKNOWN", "OTHER"),
               secondFH = rv$secondFH,
               secondFHLayer = rv$secondFHLayer,
               baseFire = rv$baseFire
    ),
    start.SEASON = rv$startTimespan,
    end.SEASON = rv$endTimespan,
    max_interval = rv$max_interval
  )
  # Save input settings to a list and then append into FH analysis object
  # FHAnalysis$AnalysisInputs<-list(
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
  # )
  print("Save input settings to a list and then append into FH analysis object")
  # FHAnalysis$FH_IDr <-
  #   fasterize(
  #     sf = FHAnalysis$OutDF,
  #     raster = cropRasters$Raster,
  #     field = "ID",
  #     fun = "first"
  #   )
  rv$FHAnalysis$FH_IDr <-
    terra::rasterize(
      x = terra::vect(rv$FHAnalysis$OutDF),
      y = eval(rv$cropRasters$rasterDef),
      field = "ID"
    )
  
  rv$cropRasters$FH_ID<-values(rv$FHAnalysis$FH_IDr)
  rv$cropRasters$FH_ID[is.na(rv$cropRasters$RGN)]<-NA
  #mask b region values ( which is masked by plm where public land only)
  print("made FHAnalysis$FH_IDr")
  
  
  
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
    rv$FHAnalysis$YSFNames <- c(FHAnalysis$YSFNames, "YSFNoBurn")
    rv$FHAnalysis$LBYNames <- c(FHAnalysis$LBYNames, "LBYNoBurn")
    rv$FHAnalysis$LFTNames <- c(FHAnalysis$LFTNames, "LFTNoBurn")
    
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
    GS_Summary <-
      makeGS_Summary(
        myFHAnalysis = rv$FHAnalysis,
        myAllCombs = rv$allCombs
      )
    
    
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
    try({
      # run relative abundance analysis  ----
      
      startBaseline <- as.integer(rv$startBaseline)
      endBaseline <- as.integer(rv$endBaseline)
      
      Baseline <- startBaseline:endBaseline
      if (rv$spListChoice == FALSE) {
        rv$TaxonList <-
          read_csv("./ReferenceTables/FAME_TAXON_LIST.csv")
      } else {
        rv$TaxonList <- read_csv(rv$customSpList)
      }
      
      if(rv$RasterRes == 75){gsub("/225m/","/75m/",rv$TaxonList$HDMPath)
      } else if(rv$RasterRes == 225){
          gsub("/225m/","/75m/",rv$TaxonList$HDMPath)}
      
      HDMPathsToCheck<-rv$TaxonList$HDMPath[rv$TaxonList$Include == "Yes"]
      MissingHDMS<-HDMPathsToCheck[!unlist(lapply(HDMPathsToCheck,file.exists))]
      if(length(MissingHDMS)>0){
        print("the following HDM rasters are missing:")
        print(MissingHDMS)
        stop("Missing HDMS")
        
          
      } 
      
      
      
      rv$HDMSpp_NO <-
        rv$TaxonList$TAXON_ID[rv$TaxonList$Include == "Yes"]
      writeSp <-
        rv$TaxonList$TAXON_ID[rv$TaxonList$WriteSpeciesRaster == "Yes"]
      writeSp <- writeSp[writeSp %in% rv$HDMSpp_NO]
      
      
      if (rv$spResponseChoice == FALSE) {
        rv$SpGSResponses <- "./ReferenceTables/OrdinalExpertLong.csv"
      } else {
        rv$SpGSResponses <-
          file.path(rv$customResponseFile)
      }
      # Select the file giving the fauna relative abundance inputs you wish to use----
      if (rv$abundByGS == TRUE) {
        AbundDataByGS <- read_csv(rv$SpGSResponses)[, c(
          "EFG_NO",
          "GS4_NO",
          "FireType",
          "Abund",
          "TAXON_ID"
        )]
        
        # If abundance data is provide by growth stage rather than time since fire expand it to the full time since fire long format ----
        rv$AbundDataLong <- AbundDataByGS %>%
          dplyr::mutate(FireTypeNo = if_else(FireType == "High", 2, if_else(FireType == "Low", 1, 0))) %>%
          dplyr::inner_join(EFG_TSF_4GS, by = c("EFG_NO", "GS4_NO")) %>%
          dplyr::arrange(TAXON_ID)
      } else {
        # Read abundance data already in full long format  ----
        rv$AbundDataLong <- read_csv(rv$SpGSResponses) %>%
          dplyr::arrange(TAXON_ID)
      }
      
      # Make the lookup list of arrays for fast calculation of cell by cell species abundance ----
      print("making Spp abund LU List")
      LU_List <- make_Spp_LU_list(
        myHDMSpp_NO = rv$HDMSpp_NO,
        myAbundDataLong = rv$AbundDataLong
      )
      
      print("finished  Spp abund LU List")
      
      print("Making spYearSumm")
      
      
      # Run the main function to get species abundance by cells ----
      rv$SpYearSumm <- calc_SpeciesRA(
        myFHAnalysis = rv$FHAnalysis,
        myAllCombs <- rv$allCombs,
        myHDMSpp_NO = rv$HDMSpp_NO,
        myWriteSpRasters = rv$makeRArasters,
        myResultsDir = rv$resultsDir,
        myLU_List = LU_List,
        myTaxonList = rv$TaxonList,
        writeYears =  rv$yearsForRasters,
        myWriteSp = writeSp,
        myCropRasters = rv$cropRasters
      )
      gc()
      # Save abundance summary outputs to csv files ----
      
      
      readr::write_csv(
        rv$SpYearSumm$SpYearSummLong,
        file.path(rv$resultsDir, "SpYearSummLong.csv")
      )
      readr::write_csv(
        rv$SpYearSumm$SpYearSummWide,
        file.path(rv$resultsDir, "SpYearSummWide.csv")
      )
      # Wrangle grouped summary species abundance depending on whether or not it has Planning unit column  ----
      
      rv$grpSpYearSummLong <- rv$SpYearSumm$grpSpYearSumm %>%
        dplyr::rename(Index_AllCombs = `myAllCombs$Index_AllCombs`) %>%
        tidyr::pivot_longer(
          -tidyr::one_of("TAXON_ID", "Index_AllCombs"),
          names_to = "SEASON",
          values_to = "sumRA"
        )
      # checks to see if planning units/burn units (PU) used
      if (match("PU", names(rv$allCombs$U_AllCombs_TFI), nomatch = 0) >
          0) {
        rv$grpSpYearSummLong <- rv$grpSpYearSummLong %>%
          dplyr::mutate(
            PU = rv$allCombs$U_AllCombs_TFI$PU[Index_AllCombs],
            EFG_NAME = rv$allCombs$U_AllCombs_TFI$EFG_NAME[Index_AllCombs]
          ) %>%
          dplyr::group_by(TAXON_ID, PU, EFG_NAME, SEASON)
      } else {
        rv$grpSpYearSummLong <- rv$grpSpYearSummLong %>%
          dplyr::mutate(EFG_NAME = rv$allCombs$U_AllCombs_TFI$EFG_NAME[Index_AllCombs]) %>%
          dplyr::group_by(TAXON_ID, EFG_NAME, SEASON)
      }
      
      rv$grpSpYearSummLong <- rv$grpSpYearSummLong %>%
        dplyr::summarise(sumRA = sum(sumRA)) %>%
        dplyr::mutate(TAXON_ID = as.integer(TAXON_ID))
      
      
      readr::write_csv(
        rv$grpSpYearSummLong,
        file.path(rv$resultsDir, "grpSpYearSummLong.csv")
      )
      
      gc()
      print("finished sp year summ")
      
      # calculate changes in abundance relative to baseline years
      rv$raDeltaAbundWide <-
        calcDeltaAbund(
          SpYearSumm = rv$SpYearSumm$SpYearSummWide,
          myFHAnalysis = rv$FHAnalysis,
          myBaseline = Baseline
        )
      readr::write_csv(
        rv$raDeltaAbundWide,
        file.path(
          rv$resultsDir,
          "SppSummChangeRelativetoBaseline.csv"
        )
      )
      gc()
      
      
      # make long form for plotting charts
      rv$raDeltaAbundLong <- rv$raDeltaAbundWide %>%
        dplyr::select(
          -tidyr::one_of(
            "TAXON_ID",
            "DIVNAME",
            "EPBC_ACT_STATUS",
            "VIC_ADVISORY_STATUS",
            "CombThreshold",
            "NoLessthanThreshhold",
            "LastLessThanThreshold"
          )
        ) %>%
        tidyr::pivot_longer(
          -tidyr::one_of(
            "COMMON_NAME",
            "SCIENTIFIC_NAME",
            "Baseline"
          ),
          names_to = "SEASON",
          values_to = "DeltaRA"
        ) %>%
        dplyr::mutate(SEASON = as.integer(SEASON))
      print("finished deltaabund")
      
      
      pathToDelete<-analysisPath
      
      analysisPath<-file.path(rv$resultsDir ,
                              paste0(rv$outputFH,
                                     ifelse(is.null(rv$customSpListName),
                                            "DefaultTaxa",
                                            gsub(".csv","",rv$customSpListName)),
                                     ".qs"))
      saveSpatRasterList(rv, analysisPath)
      unlink(pathToDelete)
      rm(analysisPath)
    })
  }
  print(paste("Saved all results for ",rv$outputFH))
  
  
  

  

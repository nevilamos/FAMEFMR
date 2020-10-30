############################################################################
# Functions used in EcoResFunctionsFMR are the calculation of Fire History
# Interval and spatial relative abundance related to fire history.
# written by nevil.amos@delwp.vic.gov.au
############################################################################


#' Generate list of species abundance lookup arrays
#' @details function creates a list of Lookup arrays for each taxon(VBA_CODE/TAXON_ID) for  YSF x EFGNO x FireTypeNo these are then used in spatial calcuation of specie abundace functions
#'
#' @param myHDMSpp_NO vector of VBA IDs for species to be incuded in analysis
#' @param myAbundDataLong long format input lookup table of species abundance x YSF xEFG_NO x FIRETYPE_NO
#'
#' @return list of 3D arrays named by TAXON_ID of relative abundance value for YSF x EFG x FIRETYPE_NO
#' @export
#'
make_Spp_LU_list <- function(myHDMSpp_NO = HDMSpp_NO,
                        myAbundDataLong = ExpertDataLong){

  myList <- list()
  for(i in myHDMSpp_NO){
    y <- myAbundDataLong[myAbundDataLong$VBA_CODE == i,]
    b = (y$YSF) + 1
    c = y$EFG_NO
    d = y$FireTypeNo
    e = y$Abund
    x <- array(NA, dim = c(max(b), 40, 4))
    for(j in 1:nrow(y)){
      x[b[j],c[j],d[j]] <- e[j]
    }
    myList[[as.character(i)]] <- x
    rm(y)
  }
  return(myList)
}

## FUNCTION makeSppYearSum2 ------------------------------------------------
# function returns summary of species summed relative abundances by year
# if writeSpRasters == TRUE, it writes species raster for each year (slows processing)
# may require excessive ram to run statewide at 75m resolution.
makeSppYearSum2 <- function(FHanalysis,
                            myHDMSpp_NO = HDMSpp_NO,
                            writeSpRasters = writeSpRasters,
                            myLU_List = LU_List,
                            ResultsDir = ResultsDir,
                            HDMVals = HDMVals,
                            TaxonList = TaxonList,
                            writeYears = NULL,
                            writeSp = writeSp) {

  # set time range for analysis
  TimeRange <- as.integer(FHanalysis$TimeSpan)
  TimeNames <- as.character(FHanalysis$TimeSpan)
  LTR <- length(TimeRange)

  # WHAT IS FH_IDr - Fire History r                                             #####-----comment-----#####
  r <- FHanalysis$FH_IDr
  # define r raster metadata
  r <- raster(nrows = nrow(r),
              ncols = ncol(r),
              ext = extent(r),
              crs = crs(r),
              vals = NULL)

  # determine what years to write out (each year or timespan)
  TimeSpan <- FHanalysis$TimeSpan
  myDF<-FHanalysis$OutDF
  if(is.null(writeYears)){
    writeYears = TimeSpan}else
    {writeYears <- writeYears[writeYears %in% TimeSpan]
    }

  # remove geometry FHanalysis DF to create a standard dataframe so columns can be subset
  # without sticky geometry of original spatial Features data frame
  st_geometry(myDF) <- NULL

  # create empty matrix with rownames containing HDM TAXON ID numbers and column names of the years (SEASONS)
  # used to house the Species abundance data for each year
  SpYearSumm <- matrix(NA, nrow = (length(myHDMSpp_NO)), ncol = LTR, dimnames = list(as.character(myHDMSpp_NO), TimeNames))

  # loop through calculation of per cell species abundance values
  # output raster values of flagged spp
  for (sp in myHDMSpp_NO) {
    cat("\r", paste( "calculating abundances for", sp))
    mySpp <- as.character(sp)

    #get the lookup array of abundance values from the list of species value lookup
    LU = myLU_List[[mySpp]]

    # gets the HDMVals for the relevant species from the matrix of HDMvalues by species.
    # values multiplied by 100 so that they can later be converted to integer if necessary without losing small values
    HDM_Vector<-as.vector(HDMVals[, grep(as.character(sp), colnames(HDMVals))]) * 100

    # makes matrices of YSF, EFG, and LFT to use in lookup from LU array
    YSF_M <- as.matrix(myDF[myAllCombs$U_AllCombs_TFI$FH_ID, FHanalysis$YSFNames]) + 1
    LFT_M <- as.matrix(myDF[myAllCombs$U_AllCombs_TFI$FH_ID, FHanalysis$LFTNames])
    EFG_M <- matrix(myAllCombs$U_AllCombs_TFI$EFG, nrow(YSF_M), ncol(YSF_M))

    # looks up the cell-wise species values for for the species abundance values by id neces in array
    Spp_M <- array(LU[cbind(as.vector(YSF_M), as.vector(EFG_M), as.vector( LFT_M))],
                   dim = dim(YSF_M))

    # Multiplies these values by cell-wise HDM values
    #effectively masking out values where the species does not occur.
    Spp_Val_Cell_Year <- Spp_M[myAllCombs$Index_AllCombs,] * HDM_Vector
    colnames(Spp_Val_Cell_Year) <- TimeNames

    # get the sum of cell values for each year for the species
    # put them in the compilation data frame
    SpYearSumm[mySpp,] <- colSums(Spp_Val_Cell_Year, na.rm = TRUE)

    #clean memory
    gc()

    # if writing species rasters is desired (function setting writeSpRasters == "Yes")
    # write out species rasters
    # this is by far the most time consuming part of the FAME processing

    if(writeSpRasters == "Yes"){
      for(myYear in as.character(writeYears)){
        cat("\r", paste("writing species abund rasters for", myYear))
        if (sp %in% writeSp|is.null(writeSp)){
          OutTif <- file.path(ResultsDir, "RA_Rasters", paste0("Spp", sp, "_", myYear, ".tif"))
          emptySpraster <- r
          values(emptySpraster) <- Spp_Val_Cell_Year[, myYear]
          writeRaster(emptySpraster,
                      OutTif,
                      options = c("COMPRESS=LZW", "TFW=YES"),
                      datatype = 'INT1U',
                      overwrite = TRUE
          )
        }
      }
    }
  } # end species raster cell abundance calcs loop

  # join species details from TaxonList  to the output tables
  SpYearSumm <- rownames_to_column(as.data.frame(SpYearSumm))
  names(SpYearSumm)[1] <- "TAXON_ID"
  TL <- TaxonList %>%
    mutate(TAXON_ID = as.character(TAXON_ID))
  SpYearSummWide <- right_join(TL, SpYearSumm)
  write.csv(SpYearSummWide,
            file.path(ResultsDir, "SpYearSummWide.csv"))
  SpYearSummLong <- right_join(TL, SpYearSumm %>%
                                 pivot_longer(-TAXON_ID, names_to = "SEASON", values_to = "SUM_RAx100"))
  write.csv(SpYearSummLong, file.path(ResultsDir, "SpYearSummLong.csv"))

  # return makeSppYearSum2 function
  return(SpYearSummWide)
}


## FUNCTION calcDeltaAbund -------------------------------------------------
# Calculates baseline RA based on input baseline years and deviation from a
# baseline (either a single year or a several years, usually a sequence
# eg 1980-1990) for each  year.
# output written to CSV "SppSummChangeRelativetoBaseline.csv"

calcDeltaAbund <- function(SpYearSummSpreadbyYear = SpYearSummWide,
                           TimeSpan = FHanalysis$TimeSpan,
                           myBaseline,
                           ResultsDir,
                           HDMSpp_NO,
                           TaxonList)
{

  # to get % of baseline need to define which columns provide the baseline
  # (one or mean of several using apply (mean)) then divide remaining values by this column.
  if (length(myBaseline == 1)){
    Baseline <- SpYearSummSpreadbyYear[, as.character(myBaseline)]
  } else {
    Baseline <- apply(SpYearSummSpreadbyYear[, as.character(myBaseline)], 1, mean)
  }

  #add baseline data to dataframe
  SpYearSummSpreadbyYear$Baseline <- Baseline

  # get integer value for current year. Used so that changes to baseline
  # are only displayed for future years or if no future years then years since baseline.
  ThisYear <- as.integer(format(Sys.Date(), "%Y"))

  SinceYear <- ifelse(sum(TimeSpan > ThisYear) > 0,
                      ThisYear,
                      max(myBaseline)
  )

  # calculate the changes from baseline
  Deltas <- as.matrix(SpYearSummSpreadbyYear[, as.character(TimeSpan[TimeSpan > SinceYear])] / Baseline)
  names(Deltas) <- paste(names(Deltas), "prop baseline")

  # calcuates two potential comparitive metrics;
  # the total number of years below threshold, and whether the last year is below threshold
  NoLessthanThreshhold <- rowSums(Deltas <= SpYearSummSpreadbyYear$CombThreshold)
  LastLessThanThreshold <- Deltas[,ncol(Deltas)] <= SpYearSummSpreadbyYear$CombThreshold

  #Subsets input columns and appends results to make an output table
  ChangeRelativeToBaseline <- cbind(SpYearSummSpreadbyYear[,c("TAXON_ID",
                                                              #"ShortName",    #####-----delelte-----######
                                                              "COMMON_NAME",
                                                              "NAME",
                                                              "DIVNAME",
                                                              "EPBC_ACT_STATUS",
                                                              "VIC_ADVISORY_STATUS",
                                                              "CombThreshold")],
                                    Deltas,
                                    NoLessthanThreshhold,
                                    LastLessThanThreshold
  )

  #writes the table to file
  write.csv(ChangeRelativeToBaseline,
            file.path(ResultsDir, "SppSummChangeRelativetoBaseline.csv"),
            row.names = FALSE)

  return(ChangeRelativeToBaseline)
}


## FUNCTION calcDraftSpList -------------------------------------------------
# Calculate the proportion of cells for the HDM in the region for each species
# works by using the indices of the standard dimesions of raster that are in
# the supplied shapefile region boundary

calcDraftSpList <- function(REG_NO,  # can this match look up table REG_LUT     #REG_NO of defined region from input (1:6) or 0 for statewide or 7 for Ad Hoc Poly),
                            RasterRes = 225,                                    ######-----'RasterRes' of 225 is defined in settings.  Change to that?
                            PUBLIC_LAND_ONLY,
                            myPoly = myPoly,       #shapefile of LF_REGIONs (default) or adhoc region,   ######-----should this be clipPoly for consistency?
                            generalRasterDir = "./InputGeneralRasters",
                            splist = "./ReferenceTables/DraftTaxonListStatewidev2.csv",
                            HDMVals = HDMVals225){                          ######-----HDMVals225 is commented out in settings, and not called anywhere prior.
  load(HDMVals)
  REG_NO <- as.integer(as.numeric(REG_NO))                                  ######-----go straight to int? rather than wrap the numeric
  splist <- read.csv(splist)
  CropDetails <- cropNAborder (REG_NO = REG_NO,
                                 RasterRes = RasterRes,
                                 PUBLIC_LAND_ONLY,
                                 myPoly = myPoly,
                                 generalRasterDir = "./InputGeneralRasters"
  )

  # get cells in polygon
  cellsInArea <- colSums(HDMVals[CropDetails$clipIDX,])
  cellsInState <- colSums(HDMVals)

  # calc proportion
  areaProp <- signif(cellsInArea / cellsInState, digits = 2)
  # pull data into dataframe
  TAXON_ID <- as.numeric(colnames(HDMVals))
  myDF <- data.frame(TAXON_ID, cellsInState, cellsInArea, areaProp)
  myDF <- left_join(splist, myDF)
  return(myDF)
}


## FUNCTION calcSppEFGLM ----------------------------------------------------
# Calculate the species in each EFG in given area for GSO calcuations.
# works by using the indices of the standard dimesions raster that are in the
# supplied shapefile region boundary (via function cropNAborder )
calcSpp_EFG_LMU <- function(REG_NO,#REG_NO of defined region from input (1:6) or 0 for statewide or 7 for Ad Hoc Poly),     ######-----sync to regions look up table?
                            RasterRes = 225,                                 ######-----'RasterRes' of 225 is defined in settings.  Change to that?
                            PUBLIC_LAND_ONLY,                        ######-----consistency for this across, T/F, Y/N
                            myPoly = myPoly,#shapefile ofLF_REGIONs( default)or  adhoc region,  ######-----should this be clippoly for consistency?
                            generalRasterDir = "./InputGeneralRasters",
                            splist = "./ReferenceTables/DraftTaxonListStatewidev2.csv",
                            HDMVals = HDMVals225,                            ######-----HDMVals225 is commented out in settings, and not called anywhere prior.
                            EFGRas = EFGRas,                                 ######-----EFGRas isn't defined anywhere prior
                            TFI_LUT = TFI_LUT){

  options(stringsAsFactors = FALSE)                                          ######-----delete? this is done at start of FAME_FMR_FHanalysisTFI_GS_calculations
  # load HDM data
  load(HDMVals)
  # read in species list
  mySpList <- read.csv(splist)[,c( "TAXON_ID","COMMON_NAME","NAME")]
  EFG <- values(raster(EFGRas))
  REG_NO <- as.integer(as.numeric(REG_NO))                                   ######-----go straight to int? rather than wrap the numeric
  CropDetails <- cropNAborder (REG_NO = REG_NO,
                                 RasterRes = RasterRes,
                                 PUBLIC_LAND_ONLY = PUBLIC_LAND_ONLY,
                                 myPoly = myPoly,
                                 generalRasterDir = "./InputGeneralRasters"
  )

  # crop EFG and HDMVals
  EFG <- EFG[CropDetails$clipIDX]
  EFG[is.na(EFG)] <- 99
  HDMVals <- HDMVals[CropDetails$clipIDX,]
  mode(EFG) <- "integer"

  # write spp EFG LMU csv
  A <- aggregate.Matrix(HDMVals, EFG ,fun = 'sum')
  myDf <- as.data.frame(as.matrix(A))
  myDf$EFG_NO <- as.integer(rownames(myDf))
  myDf <- gather(myDf, key = "TAXON_ID", "CellCount", -EFG_NO)
  myDf <- myDf[myDf$CellCount > 0,]
  myDf$TAXON_ID <- as.integer(myDf$TAXON_ID)
  myDf$EFG <- as.integer(myDf$EFG)                                              ######-----possible to bundle up the interger transforms into single line
  myDf$ha <- myDf$CellCount * 5.0625                                            ######-----cellsToHectares??  where does 5.0625 come from? cells are 225 (or 75) and neither multiply to 1ha. 5.0625 = 2.25^2 is it interchangeable if raters are 225 or 75?  --- should this be cellsToHectares function?
  myDf <- left_join(myDf, TFI_LUT[,c("EFG","EFG_NAME")], by = c("EFG_NO" = "EFG"))
  myDf <- left_join(myDf, mySpList)
  write.csv(myDf, file.path(ResultsDir, "Spp_EFG_LMU.csv"), row.names = FALSE)
  # write EFG areas csv
  EFG_AREAS <- as.data.frame(table(EFG))
  EFG_AREAS$ha <- EFG_AREAS$Freq*((RasterRes / 100) ^ 2)
  EFG_AREAS$EFG <- as.numeric(levels(EFG_AREAS$EFG))
  EFG_AREAS <- right_join(TFI_LUT[,c("EFG", "EFG_NAME")], EFG_AREAS, by = "EFG")
  write.csv(EFG_AREAS, file.path(ResultsDir, "EFG_AREAS.csv"), row.names = FALSE)
}


## FUNCTION inputRasters ----------------------------------------------------
# get names of input rasters depending on cellSize

inputRasters <- function(x = RasterRes){
  #General Input Rasters change name depending on Raster Res
  if (x == 225){
    y <- list(REGION.tif = "LF_REGION_225.tif",
              EFG.tif = "EFG_NUM_225.tif",
              PLM_GEN.tif = "PLM_GEN_225.tif",
              IDX.tif = "IndexVals225.tif",
              FIREFMZ.tif = "FIRE_FMZ_225.tif",
              DELWP.tif = "DELWP_REGION_225.tif"
    )

  }else{
    y <- list(REGION.tif = "LF_REGION_75.tif",
              EFG.tif = "EFG_NUM_75.tif",
              PLM_GEN.tif = "PLM_GEN_75.tif",
              IDX.tif = "IndexVals75.tif",
              FIREFMZ.tif = "FIRE_FMZ_75.tif",
              DELWP.tif = "DELWP_REGION_75.tif"
    )
  }
  return(y)
}


## FUNCTION makeHDMValsfromRasters  ------------------------------------------------------
# Function makes a matrix of HDM values(1,NA) constrained to those cells that
# are indexed in the cropped area
makeHDMValsfromRasters <- function(myHDMSpp_NO = HDMSpp_NO,
                                   myCropDetails = cropRasters)
{
  HDMPaths <- dir(myCropDetails$HDM_RASTER_PATH, full.names=TRUE, pattern = ".tif$")
  HDMPaths <- HDMPaths[get_Spp_No(HDMPaths) %in%
                         myHDMSpp_NO]

  print("reading HDMvalues")

  #Ncores<-8                                                                #####-----delete.  set in settings.
  # set up parallel processing
  cl<-makeCluster(Ncores)
  registerDoParallel(cl, cores=Ncores)

  myHDMVals <- foreach(i = iter(HDMPaths),.combine = cbind,.packages = "raster") %dopar% {
    myVals <- values(raster(i))[myCropDetails$IDX]
    myVals}
  stopCluster(cl)

  colnames(myHDMVals) <- as.character(get_Spp_No(HDMPaths))
  return(myHDMVals)
}

makeHDMVals <- function(myHDMSpp_NO = HDMSpp_NO,
                        myCropDetails = cropRasters,
                        RasterRes = FHanalysis$RasterRes){
  load(paste0("./HDMS/HDMVals", RasterRes, ".rdata"))
  myHDMVals <- HDMVals[myCropDetails$IDX, as.character(myHDMSpp_NO)]
  return(myHDMVals)
}

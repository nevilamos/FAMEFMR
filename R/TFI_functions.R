############################################################################
# Functions used in TFI_functionsFMRv2 are ...                                               #####-----brief description of overall function
# written by nevil.amos@delwp.vic.gov.au
############################################################################


## FUNCTION LBY_f -----------------------------------------------------------
# Function to calculate last burnt year (LBY) from matrix of rows of fire season
# iterating by year (y) used in calc_TFI_2
#' Calculate last burned year matrix (LBY)
#' @details Function to calculate last burnt year (LBY) from matrix of rows of fire season
#'  iterating by year (y) used in calc_TFI_2
#' @param M numeric matrix of fire sequences sequence in rows, values are SEASON
#' @param y numeric  SEASON
#'
#' @return matrix  of last burned year row for each unique fire history column for each SEASON
#' @import raster
#' @export
LBY_f <- function(M, y){
  M[M > y | M == 0] <- NA
  LBY <- Rfast::rowMaxs(M, value = TRUE)
  LBY[is.infinite(LBY)] <- NA
  return(LBY)
}


#' Main Tolerable fire interval (TFI) status calculation
#' @details Calculates where each cell is currently at below MinTFI or above
#'   MAX_TFI returns the per cell and long table summarised by multiple admin
#'   units and evc
#'
#' @param FHanalysis 	list containing all the fire history spatial attributes
#'   created by function fhProcess()
#' @param U_AllCombs_TFI data.table giving all combinations of cell values from
#'   the input rasters for the FAME analysis made in function calcU_All_Combs()
#' @param Index_AllCombs integer index mapping U_AllCombs_TFI to raster cells
#'   made in function calcU_All_Combs()
#' @param TFI_LUT data.frame Lookup table from EFG for
#'   "MIN_LO_TFI","MIN_HI_TFI","MAX_TFI","EFG_NAME", read form settings
#' @param OutputRasters logical whether to uptut rasters of TFI status for each
#'   year
#'
#' @return
#' @export

calc_TFI_2 <- function(FHanalysis,
                       U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI,
                       Index_AllCombs = myAllCombs$Index_AllCombs,
                       TFI_LUT,
                       OutputRasters = makeTFIRasters){


  TimeRange <- as.integer(FHanalysis$TimeSpan)
  TimeNames <- as.character(FHanalysis$TimeSpan)
  LTR <- length(TimeRange)
  r <- FHanalysis$FH_IDr #raster with spatial attributes for analysis previously stored
                          #in analysis object used as template for output rasters
  FH_ID <- values(r)
  r <- raster::raster(nrows = nrow(r),
              ncols = ncol(r),
              ext = extent(r),
              crs = crs(r),
              vals = NULL
  )
  # clean memory
  gc()

  # get the sf polygon data frame contain all the FHanalysis  attributes created in fhProcess() function
  OutTab <- FHanalysis$OutDF
  sf::st_geometry(OutTab) <- NULL
  # get ID, seasons and season types information for matrix
  ID <- OutTab$ID
  #INTFields <- names(OutTab)[grep("^INT", names(OutTab))]                      #####-----remove-----#####
  SEASFields <- names(OutTab)[grep("^SEAS", names(OutTab))]
  TYPEFields <- names(OutTab)[grep("^FireType", names(OutTab))]

  # make season/fire type matrices
  SEAS <- as.matrix(OutTab[, SEASFields]) #The Season
  SEAS[is.na(SEAS)] <- 0
  TYPE <- as.matrix(OutTab[, TYPEFields[]]) #The type of the fire

  TYPE_HI <- TYPE == 2                                                          #####-----I assume this is bushfire (type 2)? call up the look up table instead?
  TYPE_HI[TYPE_HI < 1] <- NA
  SEAS_HI <- SEAS * TYPE_HI

  TYPE_LO <- TYPE == 1                                                            #####-----I assume this is planned burn (type 1)? call up the look up table instead?
  TYPE_LO[TYPE_LO < 1 ] <- NA
  SEAS_LO <- SEAS * TYPE_LO
  # clean variable no longer needed
  rm(OutTab)
  # clean memory
  gc()

  # Calc Last Burned Year (LBY) for each firetype
  # (currently only deals with two firetypes)
  # make Lo & Hi empty dataframes for populating
  LBY_HI <- matrix(NA, nrow(SEAS), LTR)
  colnames(LBY_HI) <- TimeNames
  LBY_HI <- LBY_LO <- cbind(LBY_HI,ID)

  # loop through time range to get hi and low last burnt year
  for(i in 1:LTR){
    try({
      y = TimeRange[i]

      LBY_HI[,i] <- LBY_f(M = SEAS_HI, y)
      LBY_LO[,i] <- LBY_f(M = SEAS, y)
      cat("\r", paste("calculating LBY for", y))

    })
  }

  # partial inflation using U_AllCombs_TFI$FH_ID
  # and TFI vaules using U_AllCombs_TFI$MIN_LO_TFI
  LBY_LO <- LBY_LO[U_AllCombs_TFI$FH_ID,]
  LBY_HI <- LBY_HI[U_AllCombs_TFI$FH_ID,]

  # Calc the TFI status
  # this is the section to check if there are unusual TFI statuses
  TFI_LO <- (t(TimeRange - t(LBY_LO[,TimeNames])) - U_AllCombs_TFI$MIN_LO_TFI) < 1     #####-----for easier traceability, eplicitly declare U_All_Combs_TFI or call it from the original myAllCombs$
  TFI_HI <- (t(TimeRange - t(LBY_HI[,TimeNames])) - U_AllCombs_TFI$MIN_HI_TFI) < 1     #####-----also a short note on these numbers would be nice. e.g. 1 as it means low fire interval? why 5?
  TFI_MAX <- (t(TimeRange - t(LBY_LO[,TimeNames])) > U_AllCombs_TFI$MAX_TFI) * 5L
  TFI_COMB <- TFI_LO                                                                   #####-----can TFI_COMB just be declared where TFI_LO is?
  TFI_COMB[TFI_HI == TRUE] <- TRUE

  TFI_VAL <- TFI_COMB + TFI_MAX


  TFI_VAL[is.na(TFI_VAL)] <- -99L



  #This next section is only run for debugging unusual TFI statuses it allows their isolation at the level of unique combination of fire history and EFG
  # Check_TFI<-cbind(TFI_VAL,U_AllCombs_TFI) %>%
  #     dplyr::select(-FIRE_REG,-FIREFMZ,-PLM,-DELWP) %>%
  #       pivot_longer(all_of(TimeNames),names_to="SEASON",values_to="TFI_VAL") %>%
  #         filter(!TFI_VAL%in%c(-99,0,1,5)) %>%
  #           group_by(EFG,MIN_LO_TFI,MIN_HI_TFI,MAX_TFI,EFG_NAME,Index,FH_ID, SEASON,TFI_VAL) %>%
  #             summarize(Cells=sum(nPixel))
  #
  # Check_TFI<-data.table::as.data.table(Check_TFI)
  # setkey(Check_TFI,"FH_ID")
  # write.csv(Check_TFI,"Check_TFI.csv")
  # OutTab<-data.table::as.data.table(OutTab)
  # Check_TFI<-OutTab[Check_TFI]
  #END of DEBUGGING CODE

  # prepare output data summary tables via dplyr wrangling
  TFI_Summary <- cbind(TFI_VAL, U_AllCombs_TFI) %>%
    tidyr::pivot_longer(all_of(TimeNames), names_to = "SEASON", values_to = "TFI_VAL") %>%
    dplyr::group_by(EFG_NAME,
                    FIRE_FMZ_NAME,
                    FIRE_FMZ_SHORT_NAME,
                    FIRE_REGION_NAME,
                    DELWP_REGION,
                    EFG,
                    FIRE_REG,
                    FIREFMZ,
                    PLM,
                    DELWP,
                    SEASON,
                    TFI_VAL) %>%
    dplyr::summarize(nCells = sum(nPixel), Hectares = sum(Hectares))

  TFI_Summary <- left_join(TFI_Summary, TFI_STATUS_LUT)

  # raster output
  #Ratify all combinations raster index allows creation and export of Tif with raster attribute table
  # that can be read by ARCGIS method from : https://gis.stackexchange.com/questions/257204/saving-geotiff-from-r

  if (OutputRasters == TRUE){
    raster::values(r) <- Index_AllCombs
    rasterDatatype <- ifelse(max(Index_AllCombs) <= 65534, 'INT2S', 'INT4S') #selects the most efficient datatype depending on the size of integers in the input
    r <- raster::ratify(r)
    colnames(TFI_VAL) <- paste0("TFI_", colnames(TFI_VAL))
    levels(r)[[1]] <- cbind(levels(r)[[1]], as.data.frame(TFI_VAL))
    levels(r)[[1]] %>%
      dplyr::rename(VALUE = ID) %>%
      dplyr::mutate(VALUE = as.integer(VALUE)) %>%
      foreign::write.dbf(.,
                         file.path(ResultsDir,
                                   "TFI_Rasters",
                                   'TFI_BY_YEAR.tif.vat.dbf'),
                         factor2char = TRUE,
                         max_nchar = 254)
    raster::writeRaster(r,
                        file.path(ResultsDir,
                                  "TFI_Rasters",
                                  "TFI_BY_YEAR.tif"),
                        datatype = rasterDatatype,
                        overwrite=TRUE,
                        options=c("COMPRESS=LZW", "TFW=YES"))

  }
  return(TFI_Summary)
}

#' Calculate area BBTFI and BBTFI rasters
#' @details Calculate summary  area burned below TFI BBTFI for each SEASON in analysis
#'   (accommodating  Hi and Lo fire intensity of first burn to determine TFI)
#'   and cumulative area BBTFI.  Also optionally outputs rasters mapping areas BBTFI
#' @param myFHAnalysis list of Fhire history analysis components
#'   created by function fhProcess()
#' @param myAllCombs list made by function calc_U_AllCombs
#' @param makeBBTFIrasters logical whether or not to export
#' rasters for BBTFI to disk
#' @param myResultsDir path of directory where results will be written usually
#'   generated  by FAME script
#' @return list containing:
#' \itemize{
#' \item the date sequence matrix for each cell of the raster
#' \item the EFG TFI Lookup for each cell of the raster
#' \item the raster resolution used.
#' \item Optionally outputs rasters of BBTFI to disk if makeBBTFIrasters==TRUE.
#' }
#' @importFrom rlang .data
#' @export
calcBBTFI_2 <- function(myFHAnalysis = FHAnalysis,
                        myAllCombs = allCombs,
                        makeBBTFIrasters = makeBBTFIrasters,
                        myResultsDir = NULL)
{
  #import the allCombs data table and vector of raster cell values
  U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI
  Index_AllCombs = myAllCombs$Index_AllCombs
  # import fire history raster and make template
  r <- myFHAnalysis$FH_IDr
  FH_ID <- raster::values(r)
  r <- raster::raster(nrows = nrow(r),
              ncols = ncol(r),
              ext = raster::extent(r),
              crs = raster::crs(r),
              vals = NULL
              )


  #read the sf polygon data.frame all the fire history spatial attributes created by function fhProcess()
  #convert to a data.frame
  OutTab <- myFHAnalysis$OutDF
  sf::st_geometry(OutTab) <- NULL

  # make fire interval and type matrices
  INTFields <- names(OutTab)[grep("^INT", names(OutTab))]            # fields from FH analysis giving inter fire intervals of 1:nth fire
  SEASFields <- names(OutTab)[grep("^SEAS", names(OutTab))]          # fields giving season of 1:nth fire
  TYPEFields <- names(OutTab)[grep("^FireType", names(OutTab))]      # fields giving type of 1:nth fire


  INTS <- as.matrix(OutTab[, INTFields])                             #The fire intervals in numbered sequence
  SEAS <- as.matrix(OutTab[, SEASFields[-1]])                        #The Season of the second fire of each interval to get the season when burning under TFI occurs
  SEAS[SEAS == 0] <- NA
  TYPE <- as.matrix(OutTab[, TYPEFields[-1 * length(TYPEFields)]])   #The type of the first fire of the interval to determine whether high or lo TFI applies
  TYPE2 <- as.matrix(OutTab[, TYPEFields[-1]])                       #The FIRETYPE of the second fire - that determines whether the fire at the date that the second burn event occured was high or low for reporting
  TYPE2[TYPE2 == 0] <- NA


  TYPE_HI <- TYPE == 2                                               # TRUE where the type of the first fire is HI(2)
  TYPE_LO <- TYPE == 1                                               # TRUE where the type of the first fire is LO(1)

  # expand matrices by FH_ID
  SEAS <- SEAS[U_AllCombs_TFI$FH_ID,]
  INTS <- INTS[U_AllCombs_TFI$FH_ID,]
  TYPE_HI <- TYPE_HI[U_AllCombs_TFI$FH_ID,]
  TYPE1 <- TYPE[U_AllCombs_TFI$FH_ID,]
  TYPE2 <- TYPE2[U_AllCombs_TFI$FH_ID,]

  # make  matrices populated only with SEAS where where the inter-fire intervals are below TFI  when burned

  #if any fire interval is less than the MIN_LO_TFI then the veg has been burnt
  #below TFI multiplication of the true or false by the season returns the
  #season where true and zero where false
  BB_LO_TFI_SEASON <- SEAS * (INTS < U_AllCombs_TFI$MIN_LO_TFI)
  #remove zeros because later going to want to calculate minimum dates
  BB_LO_TFI_SEASON[BB_LO_TFI_SEASON == 0] <- NA
  #only in cases where the first FIRETYPE is HI and the interval is below
  #MIN_HI_TFI is the veg burnt below TFI
  BB_HI_TFI_SEASON <- SEAS * TYPE_HI * (INTS < U_AllCombs_TFI$MIN_HI_TFI)
  #remove zeros because later going to want to calculate minimum dates
  BB_HI_TFI_SEASON[BB_HI_TFI_SEASON == 0] <- NA

  # next two line combine dates for fires BBTFI that are below the low TFI threshold,
  # and those below the Hi TFI threshold (as determined by the TYPE of the first fire)
  BBTFI_COMB <- BB_LO_TFI_SEASON
  BBTFI_COMB[is.na(BB_LO_TFI_SEASON)] <- BB_HI_TFI_SEASON[is.na(BB_LO_TFI_SEASON)]

  # cumulative number of times BBTFI events in each combination of EFG with fire history and grouping polygons
  cumBBTFI <- t(apply(!is.na(BBTFI_COMB), 1, FUN = "cumsum"))
  totalTimesBBTFI <-rowSums(!is.na(BBTFI_COMB))

  # minimum SEASON in row determines first time that was burned below TFI
  FirstBBTFI <- apply(BBTFI_COMB, 1, min, na.rm = TRUE )
  FirstBBTFI[is.infinite(FirstBBTFI)] <- NA
  # need to set cells NA where NA in BBTFI
  cumBBTFI[is.na(BBTFI_COMB)] <- NA
  colnames(cumBBTFI) <- gsub("SEAS", "TBTFI", colnames(cumBBTFI))

  # this needs a little more annotation about what/why e.g. make 3 dataframes (WIDEish, LONG, WIDE) and their differences   #####-----are the three outputs fundamentally different? does each user use each one, or only one of them?
  # have called this combined item BBTFI_WIDEish because is is the same
  # number of rows as the true wide item but cumBBTFI is in sequence for still not spread by year.
  #BBTFI_WIDEish<-cbind(U_AllCombs_TFI,BBTFI_COMB,cumBBTFI,TYPE2,FirstBBTFI,totalTimesBBTFI)      #####-----remove
  BBTFI_WIDEish <- cbind(U_AllCombs_TFI, SEAS, cumBBTFI, TYPE2, FirstBBTFI, totalTimesBBTFI)
  BBTFI_LONG <- BBTFI_WIDEish %>%
    #select(-FirstBBTFI,-totalTimesBBTFI) %>%
                  tidyr::pivot_longer(tidyselect::all_of(names(BBTFI_WIDEish)[(ncol(U_AllCombs_TFI)+1):(ncol(BBTFI_WIDEish)-2)]),
                               names_to = c(".value", "SEQ"),
                               names_pattern = "([^0-9]+)([0-9]+)")

  BBTFI_WIDE <- BBTFI_LONG %>%
    dplyr::select(-SEQ,-totalTimesBBTFI,-FirstBBTFI) %>%
    tidyr::pivot_wider(names_from = SEAS,
                  values_from = c(TBTFI, FireType),
                  values_fn = sum,
                  names_sort = TRUE)
  #needed to deal with pivot_wider above returning list of list for non unique cases
  #using helper function unlistPivot_wider
  for (i in grep("Hectares", names(BBTFI_WIDE)):ncol(BBTFI_WIDE)){
    BBTFI_WIDE[,i] <- unlistPivot_wider(BBTFI_WIDE[,i])
  }

  BBTFI_LONG_Summary <- BBTFI_LONG %>%
    dplyr::group_by(PLM, FIRE_REGION_NAME, DELWP_REGION, EFG_NAME, FIRE_FMZ_NAME, SEAS, FireType, TBTFI) %>%
    dplyr::summarize(Hectares = sum(Hectares))

  # output rasters if applicable
  if(makeBBTFIrasters){
    print(file.path(myResultsDir,
                    "BBTFI_Rasters"))
    dir.create(file.path(myResultsDir,
                         "BBTFI_Rasters"), showWarnings = TRUE)

    raster::values(r) <- Index_AllCombs
    rasterDatatype <- ifelse(max(Index_AllCombs) <= 65534, 'INT2S', 'INT4S')
    r <- raster::ratify(r)

    levels(r)[[1]] <- cbind(levels(r)[[1]], as.data.frame(BBTFI_WIDE), FirstBBTFI, totalTimesBBTFI)
    levels(r)[[1]] %>%
      dplyr::rename(VALUE =ID) %>%
      dplyr::mutate(VALUE = as.integer(VALUE)) %>%
          foreign::write.dbf(.,
                             file.path(myResultsDir,
                                       "BBTFI_Rasters",
                                       "BBTFI_BY_YEAR.tif.vat.dbf"),
                             factor2char = TRUE,
                             max_nchar = 254)
                print(paste("BbtfiTiff=",file.path(myResultsDir, "BBTFI_Rasters", 'BBTFI_BY_YEAR.tif')))#debug
    raster::writeRaster(r,
                file.path(myResultsDir, "BBTFI_Rasters", 'BBTFI_BY_YEAR.tif'),
                datatype = rasterDatatype,
                overwrite=TRUE,
                options=c("COMPRESS=LZW", "TFW=YES")
                )


  }

  return(list("BBTFI_WIDE" = BBTFI_WIDE, "BBTFI_LONG" = BBTFI_LONG_Summary, "BBTFI_WIDEish" = BBTFI_WIDEish))

}





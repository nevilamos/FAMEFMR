#' Make long format Growth Stage Lookup matrix
#' @details expands a growth stage lookup table (provided in settings file) from four growth stages (1:4)
#'   per EFG with their years since fire spans as min(YSF) and max(YSF) to an array
#'   with YSF as row, EFG_NO as column and growth stage (1:4) as value.   NOTE:
#'   YSF has 1 added to both the Lookup and the input to deal with YSF==0 which
#'   cannot be used in the array indexing
#' @param EFG_TSF_4GS data.fame of growth stages for each EFG with start and end years
#'
#' @return matrix rows YSF, columns EFG_NO, values GS number (1:4)
#' @export
#'
makeGS_LU <- function(EFG_TSF_4GS = myEFG_TSF_4GS){
  y <- EFG_TSF_4GS
  b = (y$YSF) + 1
  c = y$EFG_NO
  e = y$GS4_NO
  x <- array(NA, dim = c(max(b), 40))
  for(j in 1:nrow(y)){
    x[b[j], c[j]] <- e[j]
  }
  return(x)
}



#' Summarise area by growth stage.
#' @param FHanalysis list containing all the fire history spatial attributes
#'   created by function fhProcess
#' @param U_AllCombs_TFI data.table giving all combinations of cell values from
#'   the input rasters for the FAME anaysis
#' @param Index_AllCombs integer index mapping U_AllCombs_TFI to raster cells #'
#' @return list of two data.frames grouped by EFG, EFG_NAME, PLM ,FIRE_FMZ_NAME,
#'   FIRE_REGION_NAME, DELWP_REGION
#'\itemize{
#'\item GS_Summary_wide Wide format table summarises area by Growth Stage and SEASON
#'\item GS_Summary_long Long format table summarises area by Growth Stage and SEASON
#'}
#' @details Generates wide and long format summary of area for each EFG and
#'   season grouped by EFG, EFG_NAME, PLM ,FIRE_FMZ_NAME, FIRE_REGION_NAME,
#'   DELWP_REGION.
#' @export
makeGS_Summary <- function(FHanalysis,
                           U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI,
                           Index_AllCombs = myAllCombs$Index_AllCombs)
  {

  TimeNames <- as.character(FHanalysis$TimeSpan)
  GS_LU <- makeGS_LU() ####----- does this change between runs? i.e. is the csv input a user input that can change.
  #####  Potentially yes so keep for now ( input .csv can be defined in settings file)

  # get the FHAnalysis sf polygon data frame ( containing all the fire history attributes) created by function fhProcess
  # and convert to a data.frame
  OutTab <- FHanalysis$OutDF
  sf::st_geometry(OutTab) <- NULL
  #ID<-OutTab$ID                                                           #####----- remove

  YSF_Fields <- names(OutTab)[grep("^YSF", names(OutTab))]
  YSFplus1 <- OutTab[, YSF_Fields] + 1
  # partial inflation using U_AllCombs_TFI$FH_ID
  YSFplus1 <- YSFplus1[U_AllCombs_TFI$FH_ID,]
  YSFplus1 <- as.matrix(YSFplus1)
  # inflate GS_LU on EFG axis using U_AllCombs_TFI$EFG
  # EFG nested within U_AllCombs_TFI in makeSppYearSum2 function
  GS_LU <- t(GS_LU)[U_AllCombs_TFI$EFG,]

  # create output summary dataframe
  indmat = cbind(rep(1:nrow(YSFplus1), ncol(YSFplus1)), as.vector(YSFplus1))
  GSrows = matrix(GS_LU[indmat], nrow = nrow(YSFplus1))
  colnames(GSrows) <- TimeNames
  GS_Summary <- cbind(U_AllCombs_TFI,GSrows)

  # populate output
  GS_Summary_Long <- GS_Summary %>%
    dplyr::select(-c(MIN_LO_TFI, MIN_HI_TFI, MAX_TFI, Index, FH_ID, FIRE_REG, FIREFMZ, DELWP)) %>%
    tidyr::pivot_longer(all_of(TimeNames),
                        names_to = "SEASON",
                        values_to = "GS") %>%
    dplyr::group_by(EFG, EFG_NAME, PLM ,FIRE_FMZ_NAME, FIRE_FMZ_SHORT_NAME, FIRE_REGION_NAME, DELWP_REGION, SEASON, GS) %>%
    dplyr::summarise(Pixels = sum(nPixel), Hectares = sum(Hectares))
  GS_Summary_Long <- left_join(GS_Summary_Long, GS_LUT)#join Growth Stage names to output table

  # change names
  setnames(GS_Summary, old = TimeNames, new = GS_Names <- paste0("GS_", TimeNames))

  return(list("GS_Summary_wide" = GS_Summary, "GS_Summary_Long" = GS_Summary_Long))
}

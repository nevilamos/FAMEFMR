#' Summarise area by growth stage.
#' @param myFHAnalysis list containing all the fire history spatial attributes
#'   created by function fhProcess
#' @param myAllCombs list made by function calc_U_AllCombs
#' @return list of two data.frames grouped by EFG, EFG_NAME, PLM ,FIRE_FMZ_NAME,
#'   FIRE_REGION_NAME, DELWP_REGION
#'\itemize{
#'\item GS_Summary_wide Wide format table summarises area by Growth Stage and SEASON
#'\item GS_Summary_long Long format table summarises area by Growth Stage and SEASON
#'}
#' @details Generates wide and long format summary of area for each EFG and
#'   season grouped by EFG, EFG_NAME, PLM ,FIRE_FMZ_NAME, FIRE_REGION_NAME,
#'   DELWP_REGION.
#' @importFrom data.table setnames
#' @export
makeGS_Summary <- function(myFHAnalysis,
                           myAllCombs)
{
  U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI
  Index_AllCombs = myAllCombs$Index_AllCombs

  TimeNames <- as.character(myFHAnalysis$TimeSpan)
  if(length(grep("NoBurn",myFHAnalysis$YSFNames))>0){TimeNames=c(TimeNames,"NoBurn")}
  GS_LU <-
    makeGS_LU() ####----- does this change between runs? i.e. is the csv input a user input that can change.
  #####  Potentially yes so keep for now ( input .csv can be defined in settings file)

  # get the FHAnalysis sf polygon data frame ( containing all the fire history attributes) created by function fhProcess
  # and convert to a data.frame
  OutTab <- myFHAnalysis$OutDF
  sf::st_geometry(OutTab) <- NULL


  YSF_Fields <- names(OutTab)[grep("^YSF", names(OutTab))]
  YSFplus1 <- OutTab[, YSF_Fields] + 1
  # partial inflation using U_AllCombs_TFI$ID
  YSFplus1 <- YSFplus1[U_AllCombs_TFI$ID, ]
  YSFplus1 <- as.matrix(YSFplus1)
  # inflate GS_LU on EFG axis using U_AllCombs_TFI$EFG
  # EFG nested within U_AllCombs_TFI in calc_SpeciesRA function
  GS_LU <- t(GS_LU)[U_AllCombs_TFI$EFG, ]

  # create output summary dataframe
  indmat = cbind(rep(1:nrow(YSFplus1), ncol(YSFplus1)), as.vector(YSFplus1))
  GSrows = matrix(GS_LU[indmat], nrow = nrow(YSFplus1))
  colnames(GSrows) <- TimeNames
  GS_Summary <- cbind(U_AllCombs_TFI, GSrows)

  # populate output
  GS_Summary_Long <- GS_Summary %>%
    dplyr::select(-c(
      MIN_LO_TFI,
      MIN_HI_TFI,
      MAX_TFI,
      Index,
      ID,
      FIRE_REG,
      FIREFMZ,
      DELWP
    )) %>%
    tidyr::pivot_longer(tidyselect::all_of(TimeNames),
                        names_to = "SEASON",
                        values_to = "GS") %>%
    dplyr::group_by(
      EFG,
      EFG_NAME,
      PLM ,
      FIRE_FMZ_NAME,
      FIRE_FMZ_SHORT_NAME,
      FIRE_REGION_NAME,
      DELWP_REGION,
      SEASON,
      GS
    ) %>%
    dplyr::summarise(Pixels = sum(nPixel), Hectares = sum(Hectares))
  GS_Summary_Long <-
    #join Growth Stage names to output table
    dplyr::left_join(GS_Summary_Long, GS_LUT)%>%
    dplyr::mutate(SEASON = as.integer(SEASON))

  # change names
  data.table::setnames(GS_Summary,
                       old = TimeNames,
                       new = GS_Names <- paste0("GS_", TimeNames))

  return(list(
    "GS_Summary_wide" = GS_Summary,
    "GS_Summary_Long" = GS_Summary_Long
  ))
}

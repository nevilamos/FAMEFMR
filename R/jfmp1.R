#' JFMP Stage 1 calculations
#'
#' @param myPUPath path to the planning unit/burn unit  shapefile covering the area.This file must contain fields:
#' \itemize{
#' \item  PU (unique integer identifier for each planning unit/burn unit)
#' \item  Hectares(area in Hectares),
#' \item  LP1_NoBurn, LP2_NoBurn, LP1_Burn, LP2_Burn - the scores for two "Life and property" metrics for each polygon in burned and unburned state at JFMPSeason0 + 4.}
#' @param grpSpYearSumm grpSpYearSumm summary of abundance of species by pivoted wide by SEASONS grouped by myAllCombs which is an output of function calc_SpeciesRA.
#' @param myAllCombs list object retuned by function calc_U_AllCombs() that contains  combinations of input raster values for the analysis
#' @param myTaxonList data.frame of species attributes (read from default or user provided .csv)
#' @param myBBTFI a list of tables produced by function calcBBTFI_2 list containing:
#' \itemize{
#' \item BBTFI_WIDE wide by SEASON  table of the number of times BBTTFI and area for each unique combination of fire history FireType,EFG, PU, and administrative subunits (District, Region etc) of area.
#' \item BBTFI_LONG long format table ( ie not spread by season) otherwise as BBTFI used for production of charts}
#' @param myJFMPSeason0 4 digit integer: the fire SEASON before the first SEASON of the JFMP being processed
#' @param zoneWt  data.frame weighting Eco vs LP for each Fire Management Zone
#' @param jfmpMetricWt data.frame weighting the four JFMP scores (Flora,Fauna,LP1 and LP2)
#' @return dataframe of planning units/ from input shapefile attributes with appended columns for Biodiversity scores ( area first BBTFI and sum of realive abundance scores weighted by the number of pixels of each species in the study area, Differences of scores between Burned and Unburned status at JFMPseason0 +4
#' @export

jfmp1 <- function(myPUPath,
                  grpSpYearSumm,
                  myAllCombs,
                  myTaxonList,
                  myBBTFI  ,
                  myJFMPSeason0,
                  zoneWt,
                  jfmpMetricWt)
{

  #Wrangle the SpYearSummRA grouped on index of all combinations of rasters, plus the TaxonList that includes count of cells in area of interest to get the weighted sum of change all species in area of interest for each PU ------

  #calcuate the SEASON after the end of the JFMP this is the SEASON for which scores are calculated
  JFMPSeason4 = myJFMPSeason0 + 4

  PU_WeightedSumRA <- grpSpYearSumm %>%
    dplyr::rename(Index_AllCombs = `myAllCombs$Index_AllCombs`) %>%
    tidyr::pivot_longer(
      -tidyr::one_of("TAXON_ID", "Index_AllCombs"),
      names_to = "SEASON",
      values_to = "sumRA"
    ) %>%
    #remove all seasons except the one for getting BURN and NO_BURN scores
    dplyr::filter(SEASON %in% c(as.character(JFMPSeason4), "NoBurn")) %>%
    dplyr::mutate(SEASON = ifelse(SEASON == "NoBurn", "NoBurn", "Burn")) %>%
    #add the PU index values
    dplyr::mutate(PU = myAllCombs$U_AllCombs_TFI$PU[Index_AllCombs]) %>%
    #group and summarise to get scores for each
    dplyr::group_by(TAXON_ID, PU, SEASON,
                    SEASON = paste0("WtSumRA_", SEASON)) %>%
    dplyr::summarise(sumRA = sum(sumRA)) %>%
    #join the cellsInArea field from taxon List
    dplyr::mutate(TAXON_ID = as.integer(TAXON_ID)) %>%
    dplyr::left_join(myTaxonList %>%
                       dplyr::select(TAXON_ID, cellsInArea)) %>%
    #calculate the for weighted sum of RA through division by
    #number of species cells (amount of species habitat) in the study area
    dplyr::mutate(weightedRA = sumRA / cellsInArea) %>%
    #sum the weighted sum of Relative abundance for all species in each PU for burned and unburned state
    dplyr::group_by(PU, SEASON) %>%
    dplyr::summarise(WeightedSumRA = sum(weightedRA)) %>%
    #calculate the difference between the sum for burned and unburned states to
    #get the score for impact of burning for the PU value will be larger
    #where burning has a negative impact. It will be negative where burning increases the abundance of a species ( smaller is good)
    tidyr::pivot_wider(names_from = SEASON, values_from = WeightedSumRA) %>%
    dplyr::mutate(WtSumRA_Diff = WtSumRA_NoBurn - WtSumRA_Burn)

  print("calculated PU weighted RA")

  #Data wrangling of BBTFI outputs to extract area BBTFI for each PU, this requires inclusion ONLY of the first time areas are burned below TFI, these are then categorised as "PastBBTFI" for all years up to JFMPSeason0, and then as BBTFI (for the first time) by the JFMP)----
  print("doing JFMPBBTFI")


  PU_BBTFI_Summ <- myBBTFI$BBTFI_LONG %>%
    # filter for first time BBTFI
    dplyr::filter(TBTFI == 1) %>%
    dplyr::mutate(JFMP_BURN = ifelse(SEAS > myJFMPSeason0, "Burn_BBTFI", "NoBurn_BBTFI")) %>%
    dplyr::group_by(PU, JFMP_BURN) %>%
    dplyr::summarise(Hectares = sum (Hectares)) %>%
    tidyr::pivot_wider(names_from = JFMP_BURN, values_from = Hectares)


  #read PU shapefile and remove geometry
  puDF <- sf::st_read(myPUPath)
  puDF$geometry <- NULL
  puDF <- puDF %>%
    #calculate the scores for LP1 and LP2 Burn minus unburned small
    #(negative) is good if LP risk is lower when burned the value will be <0
    dplyr::mutate(LP1_Diff = LP1_Burn - LP1_NoBurn,
                  LP2_Diff = LP2_Burn - LP2_NoBurn) %>%
    dplyr::left_join(PU_BBTFI_Summ) %>%
    dplyr::mutate(BBTFI_Diff = ifelse(is.na(Burn_BBTFI),0,Burn_BBTFI))%>%
    dplyr::left_join(PU_WeightedSumRA) %>%

    dplyr::mutate(WtSumRA_Diff = WtSumRA_NoBurn - WtSumRA_Burn)

  # ranking process
  puDF<-puDF%>%
    #By WHOLE AREA
    dplyr::mutate(WtSumRA_DiffStd=scales::rescale(WtSumRA_Diff,na.rm=T))%>% #standardised the Diffs ( scale 0-1) so that it makes sense to add them
    dplyr::mutate(BBTFI_DiffStd=scales::rescale(BBTFI_Diff,na.rm=T))%>% #standardised the Diffs ( scale 0-1) so that it makes sense to add them
    dplyr::mutate(LP1Std=scales::rescale(LP1_Diff,na.rm=T))%>% #standardised the Diffs ( scale 0-1) so that it makes sense to add them
    dplyr::mutate(LP2Std=scales::rescale(LP2_Diff,na.rm=T))%>% #standardised the Diffs ( scale 0-1) so that it makes sense to add them
    dplyr::mutate(WtSumRA_RankREGION=dplyr::dense_rank(WtSumRA_Diff))%>%#ranking by Diff
    dplyr::mutate(BBTFI_RankREGION=dplyr::dense_rank(BBTFI_Diff))%>%#ranking by Diff
    dplyr::mutate(LP1_RankREGION=dplyr::dense_rank(LP1_Diff))%>%#ranking by Diff
    dplyr::mutate(LP2_RankREGION=dplyr::dense_rank(LP2_Diff))%>%#ranking by Diff
    #By DISTRICT
    dplyr::group_by(DISTRICT_N)%>%
    dplyr::mutate(WtSumRA_RankDISTRICT=dplyr::dense_rank(WtSumRA_Diff))%>%#ranking by Diff within District
    dplyr::mutate(BBTFI_RankDISTRICT=dplyr::dense_rank(BBTFI_Diff))%>%#ranking by Diff within District
    dplyr::mutate(LP1_RankDISTRICT=dplyr::dense_rank(LP1_Diff))%>%#ranking by Diff within District
    dplyr::mutate(LP2_RankDISTRICT=dplyr::dense_rank(LP2_Diff))%>%#ranking by Diff within District
    #By DISTRICT and ZONE
    dplyr::group_by(DISTRICT_N,FMZ_CODE)%>%
    dplyr::mutate(WtSumRA_RankDISTRICT_FMZ=dplyr::dense_rank(WtSumRA_Diff))%>%#ranking by Diff within District and FMZ
    dplyr::mutate(BBTFI_RankDISTRICT_FMZ=dplyr::dense_rank(BBTFI_Diff))%>%#ranking by Diff within District and FMZ
    dplyr::mutate(LP1_RankDISTRICT_FMZ=dplyr::dense_rank(LP1_Diff))%>%#ranking by Diff within District and FMZ
    dplyr::mutate(LP2_RankDISTRICT_FMZ=dplyr::dense_rank(LP2_Diff))#ranking by Diff within District and FMZ

  #Scale the scores and then calculate weighted aggregate score for JFMP
  puDF<-puDF%>%dplyr::left_join(merge(zoneWt,jfmpMetricWt))%>% # calculate scaled scores with weightings to get overall score
    dplyr::mutate(DiffSum = (LP1Std * LP1Wt + LP2Std * LP2Wt) * LPwt + (WtSumRA_DiffStd * FaunaWt + BBTFI_DiffStd *FloraWt) * BDwt)

  print("finished jfmp1")

  return (puDF)
  }

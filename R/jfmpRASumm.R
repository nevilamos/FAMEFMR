
#' Summary of relative abundance changes comparing JFMP results
#'
#' @param myDraftJfmpOut data.frame of autoJFMP joined to alternative draft JFMPs mad by function joinDraftJFMP()
#' @param myGrpSpYearSummLong data.frame Long format relative abundances of species per PU and EFG and SEASON, created when calc_SpeciesRA() is run on FHAnalysis containing PU
#' @param myTaxonList data.frame of species attributes (read from default or
#'   user provided .csv)
#' @param myStartBaseline first SEASON of years to use for calculation RA baseline
#' @param myEndBaseline last SEASON of years to use for calculation RA baseline
#' @param myJFMPSeason0 the season before the commencement of the JFMP set in settings or shinyapp
#'
#' @return data.frame Summary per species for each JFMP and NoJFMP of change relative to baseline, thresholds and summed relative abundance
#' @export
jfmpRASumm <- function(myDraftJfmpOut,
                       myGrpSpYearSummLong,
                       myTaxonList,
                       myStartBaseline,
                       myEndBaseline,
                       myJFMPSeason0)
  {

  BaseLine = as.character(myStartBaseline:myEndBaseline)
  myTaxonList = myTaxonList %>%
    dplyr::select(TAXON_ID,COMMON_NAME,CombThreshold)
  jfmpNames<-c(names(myDraftJfmpOut)[grep("AutoJFMP_State",names(myDraftJfmpOut)):ncol(myDraftJfmpOut)])


  #arrange burn and NoBurn for each PU and JFMP in long format
  DF<-myDraftJfmpOut%>%
    dplyr::ungroup() %>%
    dplyr::select(tidyselect::all_of(c("PU",jfmpNames))) %>%
    dplyr::mutate(No_JFMP = ifelse(is.na(AutoJFMP_State),NA,"NO BURN")) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(c(jfmpNames,"No_JFMP")),names_to = "JFMP_Name",values_to = "Burn_NoBurn")

#calculate the baseline
  BaselineVals<-myGrpSpYearSummLong %>%
    dplyr::filter (SEASON %in% BaseLine) %>%
    dplyr::group_by(TAXON_ID,SEASON) %>%
    dplyr::summarise(sumRA = sum(sumRA)) %>%
    dplyr::group_by(TAXON_ID) %>%
    dplyr::summarise(BaselineVal = mean(sumRA))

  DF_JFMP<-myGrpSpYearSummLong %>%
    dplyr::ungroup() %>%
    dplyr::filter(SEASON %in% c(as.character(myJFMPSeason0+4),
                                "NoBurn" )) %>%
    dplyr::group_by(TAXON_ID,PU,SEASON) %>%
    dplyr::summarise(sumRA = sum(sumRA,na.rm=T)) %>%
    dplyr::mutate(Burn_NoBurn = ifelse(SEASON == "NoBurn","NO BURN","BURN")) %>%
    dplyr::select(-SEASON)


  DF2<-dplyr::left_join(DF,DF_JFMP) %>%
    dplyr::group_by(TAXON_ID,JFMP_Name) %>%
    dplyr::summarise(totalRA=sum(sumRA,na.rm=T)) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(BaselineVals) %>%
    dplyr::mutate(Delta = totalRA/BaselineVal)


  jfmpRASumm<-dplyr::right_join(myTaxonList,DF2)%>%
    dplyr::mutate(BelowThreshold = Delta<CombThreshold )

  return(jfmpRASumm)
}



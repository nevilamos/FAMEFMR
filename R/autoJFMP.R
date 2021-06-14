#' Create Automatic JFMP burn allocation
#'
#' @param dataframe of planning uit attributes and JFMP scores for each polygon burned and no t burned between JFMPSEASON0 and JFMPSEASON 0 +4 as produced by function jfmp1()
#' @param myHaTargetDF
#'
#' @return JFMP auto-nominations	 table 1, with one row for each planning unit, and columns for:
#'– Planning unit ID (x1 column)
#'– FMZ_ID category (x1)
#'– District (x1)
#'– Planning unit size in hectares (x1)
#'above columns are sourced from the input planning unit shapefile.
#'– Score (now Diff) for each metric in burn/non-burn states (2 x 4)
#'– Difference in each metric between burn/non-burn states (x4)
#'– Ranking on difference in each metric between states within District (x4)
#'– Ranking as above but within District X FMZ_ID combination (x 4)
#' @export

autoJFMP<-function(myJFMP1 = rv$puDF,myHaTargetDF = rv$JFMP_Ha_DF)
  {
  AutomJFMP_DF<-myJFMP1 %>%
    dplyr::left_join(JFMP_Ha_DF)%>%
    dplyr::arrange(DISTRICT_N,FMZ_CODE,DiffSum)%>%
    dplyr::group_by(DISTRICT_N,FMZ_CODE)%>%
    dplyr::mutate(Dist_FMZ_cum_ha=cumsum(PuAreaHa))%>%
    dplyr::mutate(cumbefore=Dist_FMZ_cum_ha-PuAreaHa)%>%
    dplyr::mutate(AutoJFMP_State=ifelse(cumbefore > TargetHa,"NO BURN","BURN"))

  return(AutomJFMP_DF)
  }

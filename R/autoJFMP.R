#' Create Automatic JFMP burn allocation
#'
#' @param  myJFMP1 data frame of planning unit attributes and JFMP scores for each polygon burned and no t burned between JFMPSEASON0 and JFMPSEASON 0 +4 as produced by function jfmp1()
#' @param myTargetHa data frame of target areas for each Fire District and FMZ_CODE in the study area in Hectares
#'
#' @return JFMP auto-nominations	 table 1, with one row for each planning unit, and columns for:
#'– Planning unit ID (x1 column)
#'– FMZ_CODE category (x1)
#'– District (x1)
#'– Planning unit size in hectares (x1)
#'above columns are sourced from the input planning unit shapefile.
#'– Score (now Diff) for each metric in burn/non-burn states (2 x 4)
#'– Difference in each metric between burn/non-burn states (x4)
#'– Ranking on difference in each metric between states within District (x4)
#'– Ranking as above but within District X FMZ_CODE combination (x 4)
#' @export

autoJFMP<-function(myJFMP1,myTargetHa)
  {
  AutomJFMP_DF<-myJFMP1 %>%
    dplyr::left_join(myTargetHa) %>%
    dplyr::arrange(DISTRICT_N,FMZ_CODE,DiffSum) %>%
    dplyr::group_by(DISTRICT_N,FMZ_CODE) %>%
    dplyr::mutate(Dist_FMZ_cum_ha=cumsum(PuHectares)) %>%
    dplyr::mutate(cumbefore=Dist_FMZ_cum_ha-PuHectares) %>%
    dplyr::mutate(AutoJFMP_State=ifelse(cumbefore > targetHa,"NO BURN","BURN"))

  return(AutomJFMP_DF)
  }

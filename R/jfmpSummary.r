#' Title Summary of JFMP areas burned and not burned
#'
#' @param myAutoJFMP JFMP table with columns for FMZ_CODE,PU,and JFMPStatus plus all the scores calculated in jfmp1()
#'
#' @return JFMP_Summary data frame CSV reporting table 2, with one row for each District and columns for:
#'– Hectares allocated to burns in auto-JFMP in each zone
#'– Total hectares allocated to burn
#'– Score for each metric (x4) if JFMP implemented
#'– Score for each metrics (x4) if JFMP not implemented
#'- Hectares allocated in each FMZ_CODE
#' @export

jfmpSummary <-
  function(myAutoJFMP = "JFMPSummary with burned or unburned") {
    JFMP_Summ <- myAutoJFMP %>%
      dplyr::ungroup() %>%
      dplyr::mutate(FaunaRAImp = if_else(AutoJFMP_State == "BURN", WtSumRA_Burn, WtSumRA_NoBurn)) %>%
      dplyr::mutate(NeverBBTFI_Imp = if_else(
        AutoJFMP_State == "BURN",
        PuHectares - (ifelse(is.na(Burn_BBTFI), 0, Burn_BBTFI)) - (ifelse(is.na(NoBurn_BBTFI), 0, NoBurn_BBTFI)),
        PuHectares - (ifelse(is.na(NoBurn_BBTFI), 0, NoBurn_BBTFI))
      )) %>%
      dplyr::mutate(LP1_Imp = if_else(AutoJFMP_State == "BURN", LP1_Burn, LP1_NoBurn)) %>%
      dplyr::mutate(LP2_Imp = if_else(AutoJFMP_State == "BURN", LP2_Burn, LP2_NoBurn)) %>%
      dplyr::group_by(DISTRICT_N) %>%
      dplyr::summarize(
        AreaHa = sum(PuHectares, na.rm = T),
        SumFaunaRAImp = sum(FaunaRAImp, na.rm = T),
        SumFaunaRANoJFMP = sum(WtSumRA_NoBurn, na.rm = T),
        SumNeverBBTFI_Imp = sum(NeverBBTFI_Imp, na.rm = T),
        SumNeverBBTFI_NoJFMP = sum(NoBurn_BBTFI, na.rm = T),
        SumLP1_Imp = sum(LP1_Imp, na.rm = T),
        SumLP1_NoJFMP = sum(LP1_NoBurn, na.rm = T),
        SumLP2_Imp = sum(LP2_Imp, na.rm = T),
        SumLP2_NoJFMP = sum(LP2_NoBurn, na.rm = T)
      )

    AreasBurnedUnburned <- myAutoJFMP %>%
      dplyr::ungroup() %>%
      dplyr::group_by(DISTRICT_N, AutoJFMP_State) %>%
      dplyr::summarise(AreaHa = sum(PuHectares, na.rm = T)) %>%
      tidyr::pivot_wider(names_from = AutoJFMP_State,
                         values_from = AreaHa) %>%
      dplyr::rename(Burned_ha = 'BURN', NoBurn_Ha = 'NO BURN')

    AreasBurnedbyFMZ_CODE <- myAutoJFMP %>%
      dplyr::ungroup() %>%
      dplyr::filter(AutoJFMP_State == 'BURN') %>%
      dplyr::group_by(DISTRICT_N, FMZ_CODE) %>%
      dplyr::summarise(AreaHa = sum(PuHectares, na.rm = T)) %>%
      tidyr::pivot_wider(names_from = FMZ_CODE,
                         values_from = AreaHa)
    names(AreasBurnedbyFMZ_CODE)[-1] <-
      paste("Burned_ha ", names(AreasBurnedbyFMZ_CODE)[-1])

    JFMP_Summ <- left_join(AreasBurnedUnburned, JFMP_Summ)
    JFMP_Summ <- left_join(JFMP_Summ, AreasBurnedbyFMZ_CODE)
    return(JFMP_Summ)

  }

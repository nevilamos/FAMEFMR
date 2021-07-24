#' Title Summary of JFMP areas burned and not burned
#'
#' @param myDraftJfmpOut JFMP table with columns for FMZ_CODE,PU,and JFMPStatus plus all the scores calculated in jfmp1()
#'
#' @return JFMP_Summary data frame CSV reporting table 2, with one row for each District and columns for:
#'– Hectares allocated to burns in auto-JFMP in each zone
#'– Total hectares allocated to burn
#'– Score for each metric (x4) if JFMP implemented
#'– Score for each metrics (x4) if JFMP not implemented
#'- Hectares allocated in each FMZ_CODE
#' @export

jfmpSummary <-
  function(myDraftJfmpOut = "JFMPSummary with burned or unburned") {


    jfmpNames<-c(names(myDraftJfmpOut)[grep("AutoJFMP_State",names(myDraftJfmpOut)):ncol(myDraftJfmpOut)])


    #arrange burn and NoBurn for each PU and JFMP in long format
    DF1<-myDraftJfmpOut%>%
      dplyr::ungroup() %>%
      dplyr::select(tidyselect::all_of(c("PU","DISTRICT_N","Burn_BBTFI",
                             "NoBurn_BBTFI","WtSumRA_Burn",
                             "WtSumRA_NoBurn","LP1_Burn",
                             "LP1_NoBurn","LP2_Burn",
                             "LP2_NoBurn","PuHectares",
                             "FMZ_CODE",jfmpNames))) %>%
      dplyr::mutate(No_JFMP = ifelse(is.na(AutoJFMP_State),NA,"NO BURN")) %>%
      dplyr::mutate(Burn_BBTFI =ifelse(is.na(Burn_BBTFI), 0, Burn_BBTFI)) %>%
      dplyr::mutate(NoBurn_BBTFI =ifelse(is.na(NoBurn_BBTFI), 0, NoBurn_BBTFI)) %>%
      tidyr::pivot_longer(cols = tidyselect::all_of(c(jfmpNames,"No_JFMP")),names_to = "JFMP_Name",values_to = "Burn_NoBurn")%>%
      dplyr::mutate(FaunaRA = ifelse(Burn_NoBurn == "BURN", WtSumRA_Burn, WtSumRA_NoBurn)) %>%
      dplyr::mutate(NeverBBTFI = ifelse(
        Burn_NoBurn == "BURN",
        PuHectares - Burn_BBTFI - NoBurn_BBTFI,
        PuHectares -  NoBurn_BBTFI)
      ) %>%
      dplyr::mutate(LP1 = ifelse(Burn_NoBurn == "BURN", LP1_Burn, LP1_NoBurn)) %>%
      dplyr::mutate(LP2 = ifelse(Burn_NoBurn == "BURN", LP2_Burn, LP2_NoBurn))

    JFMP_Summ<- DF1 %>%
      dplyr::group_by(JFMP_Name,DISTRICT_N)%>%
      dplyr::summarise(
        AreaHa = sum(PuHectares, na.rm = T),
        SumFaunaRA = sum(FaunaRA, na.rm = T),
        SumNeverBBTFI = sum(NeverBBTFI, na.rm = T),
        SumLP1 = sum(LP1, na.rm = T),
        SumLP2 = sum(LP2, na.rm = T)
      )


    AreasBurnedUnburned <- DF1 %>%
      dplyr::ungroup() %>%
      dplyr::group_by(JFMP_Name,DISTRICT_N, Burn_NoBurn) %>%
      dplyr::summarise(AreaHa = sum(PuHectares, na.rm = T)) %>%
      tidyr::pivot_wider(names_from = Burn_NoBurn,
                         values_from = AreaHa) %>%
      dplyr::rename(Burned_ha = 'BURN', NoBurn_Ha = 'NO BURN')

    AreasBurnedbyFMZ_CODE <- DF1 %>%
      dplyr::ungroup() %>%
      dplyr::filter(Burn_NoBurn == 'BURN') %>%
      dplyr::group_by(JFMP_Name,DISTRICT_N, FMZ_CODE) %>%
      dplyr::summarise(AreaHa = sum(PuHectares, na.rm = T)) %>%
      tidyr::pivot_wider(names_from = FMZ_CODE,
                         values_from = AreaHa)
    names(AreasBurnedbyFMZ_CODE)[c(-1,-2)] <-
      paste("Burned_ha ", names(AreasBurnedbyFMZ_CODE)[c(-1,-2)])

    JFMP_Summ <- dplyr::left_join(AreasBurnedUnburned, JFMP_Summ)
    JFMP_Summ <- dplyr::left_join(JFMP_Summ, AreasBurnedbyFMZ_CODE)
    return(JFMP_Summ)

  }

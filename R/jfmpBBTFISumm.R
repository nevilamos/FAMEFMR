#' Summarise areas never BBTFI by scenario.
#'
#' @param mydraftJfmpOut data.frame combining autoJFMP table and one
#' or more alternative JFMPs produced by joinDraftJFMP()
#'
#' @return data.frame summary of area never BBTFI for No Burning, AutoJFMP and each included
#'   JFMP scenario, Total for whole area and broken down by DISTRICT_N
#' @export
jfmpBBTFISumm <- function(mydraftJfmpOut) {

  jfmpNames<-c(names(mydraftJfmpOut)[grep("AutoJFMP_State",names(mydraftJfmpOut)):ncol(mydraftJfmpOut)],"No_JFMP")

  DF<-
    mydraftJfmpOut%>%
    dplyr::ungroup() %>%
    dplyr::mutate(No_JFMP = ifelse(is.na(AutoJFMP_State),NA,"NO BURN"))%>%
    tidyr::pivot_longer(cols = jfmpNames,names_to = "JFMP_Name",values_to = "Burn_NoBurn") %>%
    dplyr::mutate(
      NeverBBTFI_JFMP = ifelse(
        Burn_NoBurn == "BURN",
        PuHectares - (ifelse(is.na(Burn_BBTFI), 0, Burn_BBTFI)) - (ifelse(is.na(NoBurn_BBTFI), 0, NoBurn_BBTFI)),
        PuHectares - (ifelse(is.na(NoBurn_BBTFI), 0, NoBurn_BBTFI))
      ))

  RegSumm<- DF %>%
    dplyr::group_by(JFMP_Name) %>%
    dplyr::summarise(JFMP_TotalHaNeverBBTFI=sum(NeverBBTFI_JFMP,na.rm=T))

  Dist_Summ<-DF %>%
    dplyr::group_by(JFMP_Name,DISTRICT_N) %>%
    dplyr::summarise(JFMP_TotalHaNeverBBTFI=sum(NeverBBTFI_JFMP,na.rm=T)) %>%
    tidyr::pivot_wider(names_from = DISTRICT_N,values_from = JFMP_TotalHaNeverBBTFI)

  jfmp_bbtfi_summ <- RegSumm %>% dplyr::left_join(Dist_Summ)
}

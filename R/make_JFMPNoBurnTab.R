# Function used in  of JFMP scores calculates columns to end of FHAnalysis$OutDF for YSF scenario where there is no burning for four years from JFMPSEASON0

#' Summary of changes in relative abundance
#' @details Function used in  of JFMP scores adds columns to end of FHAnalysis$OutDF for YSF scenario wehre there is no burning for four years from JFMPSseason0
#' @param myFHAnalysis list containing all the fire history spatial attributes created by function fhProcess
#' @param JFMPSeason0 integer single SEASON (default to current SEASON) which is the season before the first year of burns ina 3 year JFMP
#' #'
#' @return data framethree columns "YSFNoBurn","LBYNoBurn","LFTNoBurn" where
#' YSFNoBurn = YSF +4 of SEASON= JFMPSeason0
#' @export
make_JFMPNoBurnTab <- function(myFHAnalysis = FHAnalysis,
                               JFMPSeason0
                               ) {
  JFMPNoBurnTab<-sf::st_set_geometry(myFHAnalysis$OutDF,value = NULL)[,paste0(c("YSF","LBY","LFT"),JFMPSeason0)]
  names(JFMPNoBurnTab)<-c("YSFNoBurn","LBYNoBurn","LFTNoBurn")
  JFMPNoBurnTab<-JFMPNoBurnTab%>%dplyr::mutate(YSFNoBurn = YSFNoBurn + 4 )
}


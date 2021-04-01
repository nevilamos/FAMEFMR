#' Calculate Geometric means of columns
#' @details Calculates geometric means of numeric columns of dataframe

#' @param x data.frame containing columns with numeric values
#' @param y vector of names of columns to exclude from calculation
#'
#' @return vector numeric or NA if not all values in column are > 0

calc_G<-function(x = raDeltaAbund, y = c("TAXON_ID", "COMMON_NAME", "SCIENTIFIC_NAME", "DIVNAME", "EPBC_ACT_STATUS",
                                         "VIC_ADVISORY_STATUS", "CombThreshold", "Baseline",  "NoLessthanThreshhold",
                                         "LastLessThanThreshold")){

  deltaabund<-x%>%select(-all_of(y))
  Gs<-apply(deltaabund,2,geoMean)
  return (Gs)
}

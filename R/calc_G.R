#' Calculate Geometric mean
#'
#' @details Calculates geometric mean of a vector of positive numeric values
#' @param x numeric vector
#'
#' @return numeric or NA if values are not all > 0
#' @export
geoMean<-function (x){
  if(any(is.na(x))){
    y = NA
    warning("vector contains NA values")
  } else if(any(x<=0)){
    y=NA
    warning("Not all x are positive values")
  } else {
    y=prod(x)^(1/length(x))
  }
  return(y)
}




#' Calculate Geometric means of columns
#'
#' @details Calculates geometric means of numeric columns of dataframe
#' @param x data.frame containing columns with numeric values
#' @param y vector of names of columns to exclude from calculation
#'
#' @return vector numeric or NA if not all values in column are > 0
#' @export
calc_G<-function(x = "raDeltaAbund", y = c("TAXON_ID", "COMMON_NAME", "SCIENTIFIC_NAME", "DIVNAME", "EPBC_ACT_STATUS",
                                         "VIC_ADVISORY_STATUS", "CombThreshold", "Baseline",  "NoLessthanThreshhold",
                                         "LastLessThanThreshold")){

  deltaabund<-x%>%dplyr::select(-tidyselect::matches(y))
  Gs<-apply(deltaabund,2,geoMean)
  Gs<-data.frame(Gs,as.integer((names(Gs))))
  names(Gs)<-c("G","SEASON")
  return (Gs)
}

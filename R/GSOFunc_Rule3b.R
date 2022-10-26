#' Rule3b: Median of survey data only
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 3b: This rule uses the median of where the observational data exists,
#' otherwise uses nothing
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns median value for the data.frame input 'Survey' column if data exists,
#' else returns NA
#' @export
#'

Rule3b <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- stats::median(datR, na.rm = TRUE)
  }
  else{
    Resp <- NA
  }
  return(Resp)
}

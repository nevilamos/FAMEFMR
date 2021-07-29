#' Rule3a: Max of survey data only
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 3a: This rule uses the maximum of where the observational data exists,
#' otherwise uses nothing
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns max value for the data.frame input 'Survey' column if data exists,
#' else returns NA
#' @export
#'

Rule3a <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- max(datR, na.rm = TRUE)
  }
  else{
    Resp <- NA
  }
  return(Resp)
}

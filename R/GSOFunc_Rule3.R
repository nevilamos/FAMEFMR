#' Rule3: Mean of survey data only
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 3: This rule uses the mean of where the observational data exists,
#' otherwise uses nothing
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns a mean value for the data.frame input 'Survey' column if data exists,
#' else returns NA
#' @export
#'

Rule3 <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- mean(datR, na.rm = TRUE)
  }
  else{
    Resp <- NA
  }
  return(Resp)
}

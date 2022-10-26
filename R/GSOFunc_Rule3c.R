#' Rule3c: Upper quartile of survey data only
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 3c: This rule uses the upper quartile of where the observational data exists,
#' otherwise uses nothing
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns upper quartile value for the data.frame input 'Survey' column if data exists,
#' else returns NA
#' @export
#'

Rule3c <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- stats::quantile(datR, probs = 0.75, na.rm = TRUE)
  }
  else{
    Resp <- NA
  }
  return(Resp)
}

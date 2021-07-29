#' Rule1: Unweighted mean survey data, otherwise mean expert data
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 1: This rule uses the mean of survey data where available, otherwise uses
#' expert opinion. It returns the mean value.
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns a mean value for the data.frame input 'Survey' column if data exists,
#' else uses 'Expert' column
#' @export
#'

Rule1 <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- mean(datR, na.rm = TRUE)
  }
  else{
    Resp <- mean(data$Response[which(data$Source == 'Expert')])
  }
  return(Resp)
}

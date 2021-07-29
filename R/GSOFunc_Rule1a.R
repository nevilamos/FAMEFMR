#' Rule1a: Max survey data where available, otherwise max expert data
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 1a: This rule uses the max of survey data where available, otherwise uses
#' expert opinion. It returns the max value.
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns a max value for the data.frame input 'Survey' column if data exists,
#' else uses 'Expert' column max
#' @export
#'

Rule1a <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- max(datR, na.rm = TRUE)
  }
  else{
    Resp <- max(data$Response[which(data$Source == 'Expert')])
  }
  return(Resp)
}

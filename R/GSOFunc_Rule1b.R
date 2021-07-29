#' Rule1b: Median survey data where available, otherwise mean expert data
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 1b: This rule uses the median of survey data where available, otherwise
#' uses expert opinion. It returns the median value.
#'
#' @param data  data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns a median value for the data.frame input 'Survey' column if data exists,
#' else uses 'Expert' column median
#' @export
#'

Rule1b <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- stats::median(datR, na.rm = TRUE)
  }
  else{
    Resp <- stats::median(data$Response[which(data$Source == 'Expert')])
  }
  return(Resp)
}

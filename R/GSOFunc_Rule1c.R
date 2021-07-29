#' Rule1c: Upper quartile survey data where available, otherwise upper quartile expert data
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 1c: This rule uses the upper quartile of survey data where available, otherwise
#' uses expert opinion. It returns the upper quartile value.
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#'
#' @return returns the upper quartile value for the data.frame input 'Survey' column if
#' data exists, else uses 'Expert' column upper quartile
#' @export
#'

Rule1c <- function(data) {
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- stats::quantile(datR, probs = 0.75, na.rm = TRUE)
  }
  else{
    Resp <-
      stats::quantile(data$Response[which(data$Source == 'Expert'), probs = 0.75])
  }
  return(Resp)
}

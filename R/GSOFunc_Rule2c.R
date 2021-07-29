#' Rule2c: Weighted average of upper quartile survey and expert data.
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' Rule 2c: This rule uses a weighted average of the upper quartile of the observational
#' data and the expert opinion where available, else uses the upper quartile of
#' the expert opinion
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#' @param Weight percentage decimal weighting towards the observational averages
#'
#' @return returns a weighted average for the upper quartile data.frame input 'Survey'
#' and 'Expert' column if 'Survey' data exists, else uses the the 'Expert'
#' column upper quartile
#' @export
#'

Rule2c <- function(data, Weight = 0.5) {
  ExpOp <-
    stats::quantile(data$Response[which(data$Source == 'Expert'), probs = 0.75])
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <-
      Weight * stats::quantile(datR, probs = 0.75, na.rm = TRUE) + (1 - Weight) * ExpOp
  }
  else{
    Resp <- ExpOp
  }
  return(Resp)
}

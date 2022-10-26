#' Rule2b: Weighted average of median survey and expert data.
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 2b: This rule uses a weighted average of the median of the observational
#' data and the expert opinion where available, else uses the median expert opinion
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#' @param Weight percentage decimal weighting towards the observational averages
#'
#' @return returns a weighted average for the median data.frame input 'Survey' and 'Expert'
#' column if 'Survey' data exists, else uses the 'Expert' column median
#' @export
#'

Rule2b <- function(data, Weight = 0.5) {
  ExpOp <- stats::median(data$Response[which(data$Source == 'Expert')])
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- Weight * stats::median(datR) + (1 - Weight) * ExpOp
  }
  else{
    Resp <- ExpOp
  }
  return(Resp)
}

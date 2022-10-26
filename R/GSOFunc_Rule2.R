#' Rule2: Weighted average of survey and expert data.
#'
#' @details Internal GSO function. Rules functions appropriately combine observational
#' and expert data. These rules are used in the functions that generates scores
#' to use in the Growth Stage Optimisation (GSO).
#' @details Rule 2: This rule uses a weighted average of the mean of the observational data
#' and the expert opinion where available, else uses the expert opinion
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#' @param Weight percentage decimal weighting towards the observational averages
#'
#' @return returns a weighted average for the mean data.frame input 'Survey' and 'Expert'
#' column if 'Survey' data exists, else uses the mean 'Expert' column
#' @export
#'

Rule2 <- function(data, Weight = 0.5) {
  ExpOp <- mean(data$Response[which(data$Source == 'Expert')])
  if (sum(data$Source == 'Survey') > 0) {
    datR <- data$Response[which(data$Source == 'Survey')]
    Resp <- Weight * mean(datR) + (1 - Weight) * ExpOp
  }
  else{
    Resp <- ExpOp
  }
  return(Resp)
}

#' Rule 0: Use mean expert opinion only
#' @details Internal GSO function. Rules functions appropriately combine
#' observational and expert data. These rules are used in the functions that
#' generates scores to use in the Growth Stage Optimisation (GSO).
#' @details Rule 0: This rule uses only expert opinion, and returns the mean value
#'
#'
#' @param data data.frame input called in OptRunCI (and subsequently used in DataGen)
#' @param Max THIS INPUT IS NOT CALLED IN FUNCTION
#'
#' @return returns a mean value for the data.frame input 'Expert' column
#' @export
#'

Rule0 <- function(data, Max = 1) {
  Resp <- mean(data$Response[which(data$Source == 'Expert')], na.rm = TRUE)
  return(Resp)
}

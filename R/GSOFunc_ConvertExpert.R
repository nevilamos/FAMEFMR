#' Convert from old to new expert opinions
#'
#' @details Non-internal GSO function.
#' Used to convert from old to new expert opinions.
#'
#'
#' @param Code TaxonID for spp. to update
#' @param x value input determining how to select updates
#' @param Update data.frame of expert scores from ./ExpertEstimate.csv
#'
#' @return returns a converted value
#' @export
#'

ConvertExpert <- function(Code, x, Update) {
  New <- Update[which(Update$TAXON_ID == Code), ]
  if (is.na(x) == TRUE) {
    y <- NA
  } else{
    if (x < 0.3) {
      y1 <- New$None
      y2 <- New$Few
      x1 <- 0
      x2 <- 0.3
    } else{
      if (x < 0.67) {
        y1 <- New$Few
        y2 <- New$Some
        x1 <- 0.3
        x2 <- 0.67
      } else{
        y1 <- New$Some
        y2 <- New$Lots
        x1 <- 0.67
        x2 <- 1
      }
    }
    m <- (y2 - y1) / (x2 - x1)
    b <- y1 - m * x1
    y <- m * x + b
  }
  return(y)
}

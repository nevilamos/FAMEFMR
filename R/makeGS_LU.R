#' Make long format Growth Stage Lookup matrix
#' @details expands a growth stage lookup table (provided in settings file) from four growth stages (1:4)
#'   per EFG with their years since fire spans as min(YSF) and max(YSF) to an array
#'   with YSF as row, EFG_NO as column and growth stage (1:4) as value.   NOTE:
#'   YSF has 1 added to both the Lookup and the input to deal with YSF==0 which
#'   cannot be used in the array indexing
#' @param myEFG_TSF_4GS data.fame of growth stages for each EFG with start and end years
#'
#' @return matrix rows YSF, columns EFG_NO, values GS number (1:4)
#' @export
#'
makeGS_LU <- function(myEFG_TSF_4GS = EFG_TSF_4GS){
  y <- myEFG_TSF_4GS
  b = (y$YSF) + 1
  c = y$EFG_NO
  e = y$GS4_NO
  x <- array(NA, dim = c(max(b), 40))
  for(j in 1:nrow(y)){
    x[b[j], c[j]] <- e[j]
  }
  return(x)
}

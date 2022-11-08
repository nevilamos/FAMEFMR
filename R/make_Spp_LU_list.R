#' Generate list of species abundance lookup arrays
#' @details function creates a list of Lookup arrays for each taxon
#'   (TAXON_ID) for  YSF x EFGNO x FireTypeNo these are then used in
#'   spatial calculation of species abundance functions
#'
#' @param myHDMSpp_NO vector of VBA IDs for species to be included in analysis
#' @param myAbundDataLong long format input lookup table of species abundance x
#'   YSF xEFG_NO x FIRETYPE_NO
#'
#' @return list of 3D arrays named by TAXON_ID of relative abundance value for
#'   YSF x EFG x FIRETYPE_NO
#' @export
make_Spp_LU_list <- function(myHDMSpp_NO = HDMSpp_NO,
                             myAbundDataLong = ExpertDataLong) {
  z <- myHDMSpp_NO %in% unique(myAbundDataLong[["TAXON_ID"]])
  if (all(z)) {
    myList <- list()
    for (i in myHDMSpp_NO) {
      y <- myAbundDataLong[myAbundDataLong$TAXON_ID == i, ]
      b = (y$YSF) + 1
      c = y$EFG_NO
      d = y$FireTypeNo
      e = y$Abund
      x <- array(NA, dim = c(max(b), 40, 4))
      for (j in 1:nrow(y)) {
        x[b[j], c[j], d[j]] <- e[j]
      }
      myList[[as.character(i)]] <- x
      rm(y)
    }
    return(myList)
  } else {
    TAXON_ID_ERROR <-
      paste(
        "Species with TAXON_IDs:",
        paste(myHDMSpp_NO[!z], collapse = " "),
        "are missing from the input Abundance data"
      )
    stop(TAXON_ID_ERROR)

  }
}

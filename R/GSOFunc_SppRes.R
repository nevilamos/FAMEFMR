#' Collate species results
#'
#' @details Internal GSO function. Used within OptRunCI()
#' @details A function to collate results for scenarios
#'
#' @param OptDat a list, the results of the DataGen() function
#' @param OptRes a data.frame, the results of the gso() function
#' @param Scen the GSO scenario (GSOScen)
#' @param TArea GSO area (GSOArea)
#'
#' @return returns a data.frame
#' @export
#'

SppRes <- function(OptDat,
                   OptRes,
                   Scen = GSOScen,
                   TArea = GSOArea) {
  res <-
    data.frame(EFG_GS = names(OptRes[-length(OptRes)]), Optimal = as.vector(unclass(t(OptRes[1, -length(OptRes)]))))
  res2 <-
    dplyr::left_join(res, Scen) %>% dplyr::left_join(TArea[, -2]) %>% dplyr::mutate(Land = PercLandscape *
                                                                                      Area) %>%
    dplyr::select(EFG_GS, Area, Scenario, Land) %>% tidyr::spread(Scenario, Land)
  res <-
    dplyr::left_join(res, res2, by = 'EFG_GS') %>% dplyr::mutate(Optimisation = Optimal *
                                                                   Area)
  return(data.frame(TAXON_ID = OptDat$Species, as.matrix(OptDat[, -1]) %*% as.matrix(res[, -1:-3])))
}

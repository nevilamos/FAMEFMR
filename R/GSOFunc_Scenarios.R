#' Joins GSO scenarios to input data
#' @details Internal GSO function.
#' @details Joins GSO scenarios to input data; prepares data to run models on scenarios
#'
#' @param data data.frame, input data
#' @param Scen data.frame, input scenarios
#' @param efg ecological fire group.  THIS INPUT IS NOT CALLED IN FUNCTION
#'
#' @return returns a data.frame containing (left) join of 'data' and 'Scen' inputs
#' @export
#'

Scenarios <- function(data, Scen = GSOScen, efg = UsedEFG) {
  Scena <- data.frame(Scenario = unique(Scen$Scenario), GMA = NA)
  for (i in 1:nrow(Scena)) {
    Scena$GMA[i] <-
      geomean.fun(dplyr::left_join(data.frame(EFG_GS = names(data)[-1]),
                                   Scen[which(Scen$Scenario == Scena$Scenario[i]), ],
                                   by = 'EFG_GS')$PercLandscape,
                  data[-1])
  }
  return(Scena)
}

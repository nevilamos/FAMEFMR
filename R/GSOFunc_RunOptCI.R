#' Run Credible Intervals on GSO calcs (main GSO function)
#'
#' @details The main Growth Stage Optimisation (GSO) function.
#' @details The Geometric Mean ABundance (GMA) for an area can be calculated using the
#' relative abundance of each species across Ecological Fire Group (EFG)
#' Growth Stages (GS) weighted by the fraction of that GS in that EFG available
#' in the area of interest. Then using non-linear optimisation, the optimal
#' allocation of GSs within each EFG can be calculated, giving a single answer,
#' based on the input from the observational data and expert opinion supplied.
#' However, there is a level of uncertainty with this data that should be
#' accounted for, as those inputs are only estimates. The uncertainty of the
#' GSO is estimated using a bootstrapping process (repeated samples and calculating 95% CIs)
#'
#' @param data data.frame, input data
#' @param efgs efgs - from wider FAME?
#' @param Scen data.frame, input scenarios
#' @param area the GSO area - wider FAME?
#' @param rule The chosen rule to apply expert opinion and observational data. Default 'Rule0'
#' @param FiTy fire type string, used as a filter in 'FireType' column of data input. Default 'High'
#' @param Weight rule weighting if a RUle 2 is to be used. Default 0.5
#' @param NRep number of repetitions of bootstrapping, with replacement. Default 20
#' @param N number of simulations
#' @param Comp The comparison species.
#'
#' @return returns a list.
#' @export
#'

OptRunCI <-
  function(data,
           efgs,
           Scen = GSOScen ,
           area = GSOArea,
           rule = 'Rule0',
           FiTy = 'High',
           Weight = 0.5,
           NRep = 20,
           N = 10,
           Comp = Comparison) {
    OptData <-
      DataGen(
        data,
        efgs = efgs,
        rule = rule,
        FT = FiTy,
        Wt = Weight,
        Rand = TRUE,
        NRep = NRep
      )[[1]]
    resultn <- gso(spp = OptData[, -1])
    resultS <- Scenarios(OptData, Scen = Scen, efg = efgs)
    resultS2 <- as.data.frame(t(c(resultS$GMA, resultn$geom)))
    names(resultS2) <-
      #c(levels(resultS$Scenario)[resultS$Scenario], 'Optimisation')
      c(resultS$Scenario, 'Optimisation')
    SpecRes <- SppRes(OptData, resultn, Scen = Scen, TArea = area)
    SpList <- array(NA, dim = c(nrow(SpecRes), ncol(SpecRes), N))
    SpList[, , 1] <- as.matrix(SpecRes)
    for (i in 2:N) {
      OptData <-
        DataGen(
          data,
          efgs = efgs,
          rule = rule,
          FT = FiTy,
          Wt = Weight,
          Rand = TRUE,
          NRep = NRep
        )[[1]]
      resultn[i, ] <- gso(spp = OptData[, -1])
      resultS <- Scenarios(OptData, Scen = Scen, efg = efgs)
      resultS2[i, ] <- t(c(resultS$GMA, resultn$geom[i]))
      SpList[, , i] <-
        as.matrix(SppRes(OptData, resultn, Scen = Scen, TArea = area))
    }
    res <- as.data.frame(t(apply(resultn, 2, CIM)))
    row.names(res)[nrow(res)] <- c('Optimisation')
    Scens <-
      as.data.frame(t(apply(resultS2 / resultS2[, which(names(resultS2) == Comp)], 2, CIM)))
    Scens0 <- as.data.frame(t(apply(resultS2, 2, CIM)))
    names(res) <- c('Prop', 'LB', 'UB')
    names(Scens) <- c('Prop', 'LB', 'UB')
    resultn2 <-
      data.frame(Scenario = rule, EFG_GS = colnames(resultn), res)[-length(resultn), ]
    resultn2$EFG <- as.numeric(stringr::str_sub(resultn2$EFG_GS, 4, 5))
    resultn2$GS <- StageNames
    resultsSp <- data.frame(apply(SpList, c(1, 2), mean))
    colnames(resultsSp) <- names(SpecRes)
    resultsSp <- dplyr::left_join(resultsSp, FaunaCodes[,c("COMMON_NAME","SCIENTIFIC_NAME","DIVNAME","TAXON_ID")])
    return(list(resultn2, res[nrow(res), ], Scens, Scens0, resultsSp))
    # return(resultn)
  }

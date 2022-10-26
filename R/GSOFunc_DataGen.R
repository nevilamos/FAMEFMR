#' Generate Growth Stage Optimisation Scores
#'
#' @details Internal GSO function. Generates the scores to use in the Growth Stage Optimisation
#' function (OptRunCI).  Includes random selection (Rand) for bootstrapping (NRep).
#'
#' @param data data.frame input called in OptRunCI
#' @param efgs efgs input called in OptRunCI
#' @param area GSO area
#' @param rule string value of internal GSO rule function to use
#' @param FT fire type string, used as a filter in 'FireType' column of data input
#' @param Wt rule weighting if a RUle 2 is to be used
#' @param Rand TRUE/FALSE to random sample for bootstrapping
#' @param NRep number of repetitions of bootstrapping, with replacement
#'
#' @return returns a list of scores
#' @export
#'

DataGen = function(data,
                   efgs,
                   area = GSOArea,
                   rule = 'Rule0',
                   FT = 'High',
                   Wt = 0.5,
                   Rand = FALSE,
                   NRep = 20) {
  nStages <- length(StageNames)
  SppList = as.data.frame(matrix(
    NA,
    nrow = length(unique(data$TAXON_ID)),
    ncol = nStages * length(efgs) + 1
  ))
  names(SppList) <-
    c('Species', paste(rep(efgs, each = nStages),
                       unique(data$GS[which(data$GS != 'NA')]),
                       sep ='_'))
  SppList$Species = unique(data$TAXON_ID)
  SppSource <- SppList
  data <- data[which(data$FireType == FT), ]
  if (Rand) {
    dataR <- data[which(data$Source == 'Expert' & data$Response >= 0), ]
    dataR$ID <- paste0('s', dataR$TAXON_ID, dataR$EFG_GS)
    dataS <-
      data[which(
        data$Source == 'Survey' &
          paste0('s', data$TAXON_ID, data$EFG_GS) %in%
          unique(dataR$ID)
      ), ]
    if (nrow(dataS) > 0) {
      dataS$ID <- paste0('s', dataS$TAXON_ID, dataS$EFG_GS)
      for (ii in 1:length(unique(dataS$EFG_GS))) {
        dataT <- dataS[which(dataS$EFG_GS == unique(dataS$EFG_GS)[ii]), ]
        Spp <- unique(dataT$TAXON_ID)
        Sites <- unique(dataT$SurvID)
        Selection <-
          sample(1:length(Sites), size = NRep, replace = TRUE)
        for (jj in 1:length(Spp)) {
          dataR <-
            rbind(dataR, dataT[which(dataT$TAXON_ID == Spp[jj]), ][Selection, ])
        }
      }
    }
    data <- dataR
  }
  for (i in 1:nrow(SppList)) {
    dat <- data[which(data$TAXON_ID == SppList$Species[i]), ]
    if (rule == 'Rule0') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule0(dat1)
        SppSource[i, j] <- 0
      }
    }
    if (rule == 'Rule1') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule1(dat1)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule1a') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule1a(dat1)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule1b') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule1b(dat1)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule1c') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule1c(dat1)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule2') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule2(dat1, Wt)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule2a') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule2a(dat1, Wt)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule2b') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule2b(dat1, Wt)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule2c') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule2c(dat1, Wt)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule3') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule3(dat1)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule3a') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule3a(dat1)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule3b') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule3b(dat1)
        SppSource[i, j] <- NA
      }
    }
    if (rule == 'Rule3c') {
      for (j in 2:ncol(SppList)) {
        dat1 <- dat[which(dat$EFG_GS == names(SppList)[j]), ]
        SppList[i, j] <- Rule3c(dat1)
        SppSource[i, j] <- NA
      }
    }
  }
  Areas <-
    dplyr::left_join(data.frame(EFG_NO = as.numeric(stringr::str_sub(
      names(SppList[, -1]), 4, 5
    )), by = "EFG_NO"), area[, -2])$Area
  SppList[, -1] <- sweep(SppList[, -1], MARGIN = 2, Areas, '*')
  NoZero <- which(rowSums(SppList[, -1], na.rm = TRUE) != 0)
  SppList[, -1] <-
    apply(SppList[, -1], c(1, 2), function(x) {
      ifelse(is.na(x), 0, x)
    })
  return(list(SppList[NoZero, ], SppSource[NoZero, ]))
}

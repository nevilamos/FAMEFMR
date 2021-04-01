# Calculates baseline RA based on input baseline years and deviation from a
# baseline (either a single year or a several years, usually a sequence
# eg 1980-1990) for each  year.
# output written to CSV "SppSummChangeRelativetoBaseline.csv"

#' Summary of changes in relative abundance
#' @details Calculates the change in relative abundance compared to a baseline SEASON or mean of SEASONS
#' @param SpYearSumm data.frame output by function calc_SpeciesRA()
#' @param myFHAnalysis list containing all the fire history spatial attributes created by function fhProcess
#' @param myBaseline integer single SEASON or sequence of SEASONS used to create the baseline relative species abundance for comparison of change
#' #'
#' @return data frame wide format summary change in relative abundance of species by SEASON relative to a Baseline
#' @export
calcDeltaAbund <- function(SpYearSumm = SpYearSummWide,
                           myFHAnalysis,
                           myBaseline
)
{
  TimeSpan = myFHAnalysis$TimeSpan
  # to get % of baseline need to define which columns provide the baseline
  # (one or mean of several using apply (mean)) then divide remaining values by this column.
  if (length(myBaseline) == 1){
    Baseline <- SpYearSumm[, as.character(myBaseline)]
    BaselineText<- as.character(myBaseline)
  } else {
    Baseline <- apply(SpYearSumm[, as.character(myBaseline)], 1, mean)
    BaselineText<- paste(min(myBaseline),"-", max(myBaseline))
  }

  #add baseline data to dataframe
  SpYearSumm$Baseline <- Baseline

  # get integer value for current year. Used so that changes to baseline
  # are only displayed for future years or if no future years then years since baseline.
  ThisYear <- as.integer(format(Sys.Date(), "%Y"))

  # SinceYear <- ifelse(sum(TimeSpan > ThisYear) > 0,
  #                     ThisYear,
  #                     max(myBaseline)
  # )
  SinceYear <- max(myBaseline)


  # calculate the changes from baseline
  Deltas <- as.matrix(SpYearSumm[, as.character(TimeSpan[TimeSpan > SinceYear])] / Baseline)
  names(Deltas) <- paste(names(Deltas), "prop baseline")

  # calculates two potential comparative metrics;
  # the total number of years below threshold, and whether the last year is below threshold
  NoLessthanThreshhold <- rowSums(Deltas <= SpYearSumm$CombThreshold)
  LastLessThanThreshold <- Deltas[,ncol(Deltas)] <= SpYearSumm$CombThreshold

  #Subsets input columns and appends results to make an output table
  SpYearSumm<-SpYearSumm[,c("TAXON_ID",
                "COMMON_NAME",
                "SCIENTIFIC_NAME",
                "DIVNAME",
                "EPBC_ACT_STATUS",
                "VIC_ADVISORY_STATUS",
                "CombThreshold")]
  SpYearSumm$Baseline <- BaselineText
  ChangeRelativeToBaseline <- cbind(SpYearSumm,
                                    Deltas,
                                    NoLessthanThreshhold,
                                    LastLessThanThreshold
  )
  return(ChangeRelativeToBaseline)
}

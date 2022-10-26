# Calculates baseline RA based on input baseline years and deviation from a
# baseline (either a single year or a several years, usually a sequence
# eg 1980-1990) for each  year.
# output written to CSV "SppSummChangeRelativetoBaseline.csv"

#' Summary of changes in relative abundance
#' @details Calculates the change in relative abundance compared to a baseline SEASON or mean of SEASONS
#' @param SpYearSumm wide format data.frame output by function calc_SpeciesRA() ( one of several
#' data.frames in the list returned)
#' @param myFHAnalysis list containing all the fire history spatial attributes created by function fhProcess
#' @param myBaseline integer single SEASON or sequence of SEASONS used to create the baseline relative species abundance for comparison of change
#' #'
#' @return data frame wide format summary change in relative abundance of species by SEASON relative to a Baseline
#' @export
calcDeltaAbund <- function(SpYearSumm = SpYearSumm$SpYearSummWide,
                           myFHAnalysis,
                           myBaseline
)
{
  TimeSpan = myFHAnalysis$TimeSpan
  # to get % of baseline need to define which columns provide the baseline
  # (one or mean of several using apply (mean)) then divide remaining values by this column.
  if (length(myBaseline) == 1){
    BaselineVals <- SpYearSumm[, as.character(myBaseline)]
    BaselineText<- as.character(myBaseline)
  } else {
    BaselineVals <- apply(SpYearSumm[, as.character(myBaseline)], 1, mean)
    BaselineText<- paste(min(myBaseline),"-", max(myBaseline))
  }

  #add baseline SEASON(s) to dataframe
  SpYearSumm$Baseline <- BaselineText

  SinceYear <- max(myBaseline)

  # calculate the changes from baseline, if clause deals with new cases where BaselineVals are tibble or data frame
  if(is.data.frame(BaselineVals)){BaselineVals <- dplyr::pull(BaselineVals)}
  Deltas <- as.matrix(SpYearSumm[, as.character(TimeSpan[TimeSpan > SinceYear])] / BaselineVals)


  # calculates two potential comparative metrics;
  # the total number of years below threshold, and whether the last year is below threshold
  NoLessthanThreshhold <- rowSums(Deltas <= SpYearSumm$CombThreshold)
  LastLessThanThreshold <- Deltas[,ncol(Deltas)] <= SpYearSumm$CombThreshold

#puts together the "labelling columns from input table with the Deltas and other results for output
  ChangeRelativeToBaseline <- cbind(SpYearSumm%>%
                                      dplyr::select(TAXON_ID,
                                                    COMMON_NAME,
                                                    SCIENTIFIC_NAME,
                                                    EPBC_ACT_STATUS,
                                                    VIC_ADVISORY_STATUS,
                                                    CombThreshold,
                                                    Baseline),
                                    Deltas,
                                    NoLessthanThreshhold,
                                    LastLessThanThreshold
  )

}

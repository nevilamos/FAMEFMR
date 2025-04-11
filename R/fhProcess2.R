#' Second stage of Fire history processing for FAME
#'
#' @description
#' Takes the results of the main FH process function  and calculates
#'   Time Since fire (TSF) and Last Fire Type (LFT) and Last burnt year (LBY)
#'   for each SEASON as defined in the input arguments, these values are appended
#'   to the output sf polygon dataframe.
#'
#'
#' @param inFH1 input partially processed Fire History from
#' @param start.SEASON first season for which output is wanted
#' (four digit year as integer), if NUll then second season in in history is used
#'  (cannot use first season because it has no interval, this may still fail if
#'   there is no overlap)
#' @param end.SEASON last season required, if NULL then maximum value in fire
#' history scenario used
#' @param max_interval maximum inter fire interval after a HIGH fire for which
#' subsequent fire is reassigned to HIGH, if 0 then no reassignment)
#'
#' @return A list containing: \itemize{
#' \item OutDF A completed Fire Analysis Simple feature collection of polygons and/or
#'  multipolygons  with columns for Fire Season sequences, inter-fire intervals
#'  Last fire type,last fire SEASON and last burned SEASON for each season selected
#'  in the input range.
#' \item TimeSpan integer vector sequence of SEASONS to include in the analysis
#'  output
#' \item YSFNames names of TSF years in output, needed by downstream functions
#' \item LBYNames names of LBY years in output, needed by downstream functions
#' \item LFTNames names of LBY years in output, needed by downstream functions }
#' @export
#'
#' @examples
#' # randomFH<-generate_random_fire_history(20)
#' # plot(randomFH)
#' # outFH1<-fhProcess1(randomFH)
#' # plot(outFH1)
#' # #complete all fields in FH analysis
#' # outFH2<-fhProcess2(outFH1)
#' #
#' # plot(outFH2$OutDF,max.plot = 20)
#'
#'

fhProcess2<-function(inFH1,
                     start.SEASON = NULL,# first season for which output is wanted (four digit year as integer), if NULL then second season in in history is used (cannot use first season because it has no interval, this may still fail if there is no overlap)
                     end.SEASON = NULL,# last season required, if NULL then largest value in fire history scenario used
                     max_interval = 0 ){# maximum inter fire interval after a HIGH fire for which subsequent fire is reassigned to HIGH, if 0 then no reassignment)
  mySF<-inFH1

  mySF<-inFH1 %>%
    dplyr::group_by(across(-geometry)) %>%
    dplyr::summarise()
  # get just the names for the season and FireType sequence columns and make sure they are in ascending order
  SEASNames <- sort(names(mySF)[grep(pattern = "SEAS", names(mySF))])
  FTNames <- sort(names(mySF)[grep(pattern = "FireType", names(mySF))])

  # calculates all inter-fire intervals note by offsetting matrices by one column
  myDF<-mySF %>%
    sf::st_drop_geometry()

  SEAS_Matrix <- as.matrix(myDF[, SEASNames])

  FT_matrix <- as.matrix(myDF[, FTNames])

  Cols <- ncol(SEAS_Matrix)

  SEAS_Matrix[SEAS_Matrix == 0] <- NA

  Interval <- SEAS_Matrix[,2:Cols] - SEAS_Matrix[,1:Cols-1]

  IntNames <- paste("INT", sprintf("%02d", 1:(Cols-1)), sep = "")
  colnames(Interval) <- IntNames

  # bind the intervals calculated back to the main dataframe
  mySF <- cbind(mySF, Interval)[,c("geometry",SEASNames,FTNames,IntNames)]

  # work around for low firetype following after high firetype
  if(max_interval>0){
    FT_matrix<-fireTypeLowToHigh(max_interval = as.integer(max_interval),Interval_Matrix = Interval,Firetype_Matrix = FT_matrix)
    mySF[, FTNames]<-FT_matrix
  }else if (max_interval<0)  {
    stop ("max interval cannot be less than 0")
  }else {

  }





  mySF <- sf::st_as_sf(mySF)

  # time span (range of consecutive years) for which Fire History sequences are calculated
  # from start and end seasons defined in the settings file
  # else if NULL; start and end season for calculation is determined by the second smallest and max values in the dataset
  # the earliest start season is the date of a second season in the input data , next 3 lines prevent manual start
  # season setting less than this value, which would cause an error.
  uniqueSEASONS<-sort(unique(as.vector(SEAS_Matrix)))
  min.SEASON<-uniqueSEASONS[2] # second season in fire history
  if(is.null(start.SEASON)){
    start.SEASON = min.SEASON
  } else {
    if(start.SEASON < min.SEASON){
      start.SEASON = min.SEASON
    } else {
      start.SEASON = start.SEASON
    }
  }

  if(is.null(end.SEASON)){
    max.SEASON <- max(uniqueSEASONS)
  } else {
    max.SEASON = end.SEASON
  }

  TimeSpan<-start.SEASON:max.SEASON



  # calculating the last burnt season for each sequences for each year
  # this process is duplicated here and in CALC_TFI2 function
  # ideally it should only be run once but it is pretty fast so probably does not matter too much.
  LTR <- length(TimeSpan)
  SEAS_Matrix[is.na(SEAS_Matrix)] <- 0
  LBY <- matrix(NA, nrow(SEAS_Matrix), LTR)
  for(i in 1:LTR){
    try({
      y = TimeSpan[i]
      LBY[,i] <- LBY_f(M = SEAS_Matrix, y)
      #print(y)
    })
  }

  # subtraction leads to matrix transposition which needs to be re transposed in the next couple of lines
  tYSF <- TimeSpan-t(LBY)
  YSF <- t(tYSF)

  # name matrix columns
  YSFNames <- paste0("YSF", TimeSpan)
  LBYNames <- paste0("LBY", TimeSpan)
  LFTNames <- paste0("LFT", TimeSpan)
  colnames(YSF) <- YSFNames
  colnames(LBY) <- LBYNames

  # Create matrix used for getting last firetype by year
  print("calculating lookup matrix for getting last FireType by SEASON")
  SEAS_Matrix[SEAS_Matrix == 0] <- NA
  LUM <- matrix(NA, nrow(SEAS_Matrix), max.SEASON)
  for (i in 1:nrow(SEAS_Matrix)){
    R <- i
    C <- as.numeric(stats::na.omit(SEAS_Matrix[i,]))
    V <- (FT_matrix[i,(1:length(C))])
    LUM[R,C] <- V
  }


  # Calculate last fire type
  print("calculating last fire type")
  LFT <- matrix(NA, nrow(SEAS_Matrix), LTR)
  for(i in 1:nrow(SEAS_Matrix)){
    LFT[i,] <- LUM[i, LBY[i,]]
  }
  colnames(LFT) <- LFTNames


  # put together output dataframe
  mySF <- cbind(mySF, YSF)
  mySF <- cbind(mySF, LBY)
  mySF <- cbind(mySF, LFT)
  mySF$ID <- as.integer(rownames(mySF))
  print("Completed making FH object")

  return(list("OutDF" = mySF,
              "YSFNames" = YSFNames,
              "LBYNames" = LBYNames,
              "LFTNames" =LFTNames,
              "TimeSpan"=TimeSpan))
}







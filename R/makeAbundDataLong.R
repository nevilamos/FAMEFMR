#' Makes Long format  Fauna  abundance table
#'
#' @details Supporting function to Make long format table for fauna abundance scores
#'  by "TAXON_ID"  ,FireType EFG and Growth Stage from input wide format table currently
#'   deals only with FireType of "High" and "Low"
#'   which are converted to 2 and 1 respectively
#' @param AbundDataByGSFile .csv input file containing fields
#' "EFG_NO", "GS4_NO",  "FireType" , "Abund", "TAXON_ID"
#' with Abund values for' both FireTypes for each growth stage "GS4_NO"
#' @param myEFG_TSF_4GS table of each combination of "EFG_NO", "GS4_NO", and "YSF" generally read in at beginning of FAME in settings file
#' @return long format table with one row for each combination of "EFG_NO", "GS4_NO",  "FireType" , "Abund", "TAXON_ID" and "YSF"
#' sorted by TAXON_ID
#' @export

makeAbundDataLong<- function(AbundDataByGSFile = "./ReferenceTables/OrdinalExpertLong.csv",
                             myEFG_TSF_4GS = EFG_TSF_4GS){
  AbundDataByGS <- utils::read.csv(AbundDataByGSFile)[,c("EFG_NO", "GS4_NO",  "FireType" , "Abund", "TAXON_ID")]
  AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="High"]<-2
  AbundDataByGS$FireTypeNo[AbundDataByGS$FireType=="Low"]<-1
  AbundDataByGS<-AbundDataByGS[!is.na(AbundDataByGS$Abund),c("EFG_NO", "GS4_NO",  "FireTypeNo" , "Abund", "TAXON_ID")]


  AbundDataLong = merge(AbundDataByGS, myEFG_TSF_4GS,   by=c('EFG_NO','GS4_NO'))
  AbundDataLong<-AbundDataLong[order(AbundDataLong$TAXON_ID),]
  return(AbundDataLong)
}

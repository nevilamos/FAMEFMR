#' Title
#'
#' @param DraftJFMPFile path to csv file containing two columns:
#'  PU  and Draft_JFMP_State which is the burn state 
#'  ("Burn" or "NoBurn") for the planning unit(PU) in the draft JFMP
#'  The PU values must be idetical to those in myAutoJFMP
#' @param AutoJFMP autoJFMP for the corresponding fire History and
#'  planning units
#'
#' @return data.frame JFMP table with one row for each planning unit,
#'  and columns for:
#' – Planning unit ID (x1 column)
#' – FMZ category (x1)
#' – District (x1)
#' – Planning unit size in hectares (x1)
#' – Burn/non-burn state from draft JFMP (x 1)
#' – Score for each metric in burn/non-burn states (2 x 4)
#' – Difference in each metric between burn/non-burn states (x4)
#' – Ranking on difference in each metric between states within District (x4)
#' – Ranking as above but within District X FMZ combination (x 4)
#'#'  This is the same as table 2.1(autoJFMP) and 2.2 except different areas are considered for burning 
#' @export

processDraftJFMP <- function(DraftJFMPFile =rv$draftJFMPFile,AutoJFMP =rv$autoJFMP) {
  # 3. JFMP Review Tool
  draftJFMP<-read_csv(DraftJFMPFile)
  
  if(!identical(sort(draftJFMP$PU),sort(AutoJFMP$PU))){
    stop("The Draft JFMP input PU \n must be identical to the\n AutoJFMP_DF PU, NA values can be assigned to those in other districts")
  }else{
    # JFMP table 1, with one row for each planning unit, and columns for:
    
    
    myJFMP_DF<-AutoJFMP %>% 
      left_join(draftJFMP) %>% 
      mutate(AutoJFMP_State = Draft_JFMP_State)
    
  }
  return(myJFMP_DF)
}
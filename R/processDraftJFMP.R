#' Process draft JFMP inputs
#' @details Takes a table of burn units with identical PU ids to the autoJFMP can makes table of burn status for the draft JFMP and relevant ecological resilience and life and property risk scores per burn unit.

#' @param DraftJFMPFile path to csv file containing two columns:
#'  PU  and Draft_JFMP_State which is the burn state
#'  ("Burn" or "NoBurn") for the planning unit(PU) in the draft JFMP
#'  The PU values must be identical to those in myAutoJFMP
#' @param AutoJFMP autoJFMP for the corresponding fire History and
#'  planning units
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
#'  This is the same as autoJFMP table  except different areas are considered for burning
#' @export
processDraftJFMP <- function(
  myDraftJFMPFile =rv$draftJFMPFile,
  myAutoJFMP =rv$autoJFMP) {

  draftJFMP<-read_csv(myDraftJFMPFile)

  if(!identical(sort(draftJFMP$PU),sort(myAutoJFMP$PU))){
    stop("The Draft JFMP input PU \n must be identical to the\n AutoJFMP_DF PU, NA values can be assigned to those in other districts")
  }else{

    myJFMP_DF<-myAutoJFMP %>%
      left_join(draftJFMP) %>%
      mutate(AutoJFMP_State = Draft_JFMP_State)

  }
  return(myJFMP_DF)
}


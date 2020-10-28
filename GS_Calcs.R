############################################################################
# Functions used in GS_Calcs are ...                                        #####-----brief description of overall function
# written by nevil.amos@delwp.vic.gov.au
############################################################################


## FUNCTION makeGS_LU ------------------------------------------------------
# Function makes LU matrix for Growth stage from TSF and EFG 
# YSF has 1 added to both the Lookup and the input to deal with YSF==0
# which cannot be used in the array indexing
makeGS_LU <- function(EFG_TSF_4GS = myEFG_TSF_4GS){
  y <- EFG_TSF_4GS
  b = (y$YSF) + 1
  c = y$EFG_NO
  e = y$GS4_NO
  x <- array(NA, dim = c(max(b), 40))
  for(j in 1:nrow(y)){
    x[b[j], c[j]] <- e[j]
  }
  return(x)
}


## FUNCTION makeGS_Sum -----------------------------------------------------
# function to sum .....                                                     #####-----add comment on what GS stands for

makeGS_Sum <- function(writeGSRasters,                                      #####----- declare the function variable, then give it the existing variable. on all
                       ResultsDir,
                       U_AllCombs_TFI = myAllCombs$U_AllCombs_TFI,
                       Index_AllCombs = myAllCombs$Index_AllCombs,
                       FHanalysis,
                       writeYears = NULL) {
  
  TimeNames <- as.character(FHanalysis$TimeSpan)
  GS_LU <- makeGS_LU()                                                     #####----- does this change between runs? i.e. is the csv input a user input that can change.  
  
  # read the corresponding shapefile for the FireHat analysis chosen
  # and convert to a data.frame
  OutTab <- FHanalysis$OutDF
  st_geometry(OutTab) <- NULL
  #ID<-OutTab$ID                                                           #####----- remove
  
  YSF_Fields <- names(OutTab)[grep("^YSF", names(OutTab))]
  YSFplus1 <- OutTab[, YSF_Fields] + 1
  # partial inflation using U_AllCombs_TFI$FH_ID
  YSFplus1 <- YSFplus1[U_AllCombs_TFI$FH_ID,]
  YSFplus1 <- as.matrix(YSFplus1)
  # inflate GS_LU on EFG axis using U_AllCombs_TFI$EFG
  # EFG nested within U_AllCombs_TFI in makeSppYearSum2 function
  GS_LU <- t(GS_LU)[U_AllCombs_TFI$EFG,]                                   

  # create output summary dataframe
  indmat = cbind(rep(1:nrow(YSFplus1), ncol(YSFplus1)), as.vector(YSFplus1))
  GSrows = matrix(GS_LU[indmat], nrow = nrow(YSFplus1))
  colnames(GSrows) <- TimeNames
  GS_Summary <- cbind(U_AllCombs_TFI,GSrows)
  
  # populate output  
  GS_Summary_Long <- GS_Summary %>%
    select(-c(MIN_LO_TFI, MIN_HI_TFI, MAX_TFI, Index, FH_ID, FIRE_REG, FIREFMZ, DELWP)) %>% # remvoe these columns from U_AllCombs_TFI
      pivot_longer(all_of(TimeNames), names_to = "SEASON", values_to = "GS") %>%
        group_by(EFG, EFG_NAME, PLM ,FIRE_FMZ_NAME, FIRE_FMZ_SHORT_NAME, FIRE_REGION_NAME, DELWP_REGION, SEASON, GS) %>%
          summarise(Pixels = sum(nPixel), Hectares = sum(Hectares))
  GS_Summary_Long <- left_join(GS_Summary_Long, GS_LUT)
  
  # change names
  setnames(GS_Summary, old = TimeNames, new = GS_Names <- paste0("GS_", TimeNames))
  
  return(list("GS_Summary_wide" = GS_Summary, "GS_Summary_Long" = GS_Summary_Long))
}
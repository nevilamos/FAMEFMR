makeGS_Sum<-function(TimeSpan = FHanalysis$TimeSpan,
                     #myHDMSpp_NO = HDMSpp_NO,
                     writeGSRasters,
                     myLU = GS_LU,
                     myResultsDir = GS_Dir,
                     EFG = cropRasters$EFG,
                     myCropDetails = cropRasters,
                     myFHResults = FHanalysis,
                     myYSF_LFT = tsf_ysf_mat,
                     writeYears=NULL) {
  GSYearSumm <- NULL
  for (year in TimeSpan) {
    print(paste("Starting GS for",year))
    myYSF <- paste0("YSF", year)
    YSF <- myYSF_LFT[, myYSF] + 1
    myDim <- length(YSF)
    Mask_idx <- (1:myDim)
    RegMaskVal <- YSF + EFG    #+RGN    
    M<-cbind(YSF,EFG)
    
    
    LU = myLU
    
    getVals <- Mask_idx[!is.na(RegMaskVal)]#
    OutTif <-file.path(myResultsDir, paste0("GS_YR_", year, ".tif"))
    #OutName<-paste0(year,"_",sp)
    
    Out <- array(NA, myDim)
    
    
    Out[getVals]<-as.integer(LU[M[getVals,]])
    
    if(writeGSRasters=="Yes"){
      
      
      if(year%in%writeYears|is.null(writeYears)){
        
        emptyraster <- myCropDetails$Raster
        values(emptyraster) <- as.vector(Out)
        writeRaster(
          emptyraster,
          OutTif,
          options = c("COMPRESS=LZW", "TFW=YES"),
          datatype = 'INT1U',
          overwrite = TRUE
        )
        
      }
    }
    GSYrres <- c(year, table(Out))
    GSYearSumm <- rbind(GSYearSumm, GSYrres)
    
    print(paste("Finished",year))
    
  }
  return(GSYearSumm)
}

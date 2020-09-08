###makeSppYearSum#######################################################
#function returns summary of species summed relative abundances by year
#also if writeSpRasters==TRUE it writes species raster sor each year (slows processing)
#works as fast as old foreach version without blowing out RAM for statewide 225m
makeSppYearSum2<-function(FHanalysis,
                         myHDMSpp_NO = HDMSpp_NO,
                         writeSpRasters = writeSpRasters,
                         myLU_List = LU_List,
                         #YSF_TSF_Dir = YSF_TSF_Dir,
                         ResultsDir = ResultsDir,
                         #EFG = cropRasters$EFG,
                         myCropDetails = cropRasters,
                         HDMVals = HDMVals,
                         myFHResults = FHanalysis,
                         #myYSF_LFT = tsf_ysf_mat,
                         TaxonList = myTaxonList,
                         writeYears=NULL,
                         writeSp =NULL,
                         ) {
  
  #create NULL object to house the Species abundance data for each year.
  TimeRange<-as.integer(FHanalysis$TimeSpan)
  TimeNames<-as.character(FHanalysis$TimeSpan)
  LTR<-length(TimeRange)
  r<-FHanalysis$FH_IDr
  
  r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
  
  
  
  
  
  SpYearSumm <- NULL
  TimeSpan<-FHanalysis$TimeSpan
  myDF<-FHanalysis$OutDF
  st_geometry(myDF)<-NULL
  
  
  SpYearSumm<-matrix(NA,nrow=(length(myHDMSpp_NO)),ncol=LTR,dimnames=list(as.character(myHDMSpp_NO),TimeNames))
  for (sp in myHDMSpp_NO) {
    tic(sp)
    mySpp<-as.character(sp)
    LU = myLU_List[[mySpp]]
    HDM_Multiplier<-HDMVals[, ]
    YSF_M<-as.matrix(myDF[myAllCombs$U_AllCombs_TFI$FH_ID,FHanalysis$YSFNames])+1
    LFT_M<-as.matrix(myDF[myAllCombs$U_AllCombs_TFI$FH_ID,FHanalysis$LFTNames])
    EFG_M<-matrix(myAllCombs$U_AllCombs_TFI$EFG,nrow(YSF_M),ncol(YSF_M))
    
    Spp_M<-array(LU[cbind(as.vector(YSF_M),as.vector(EFG_M),as.vector( LFT_M))],dim=dim(YSF_M))
    Spp_Val_Cell_Year<-Spp_M[myAllCombs$Index_AllCombs,]*HDMVals[, as.character(sp)]
    colnames(Spp_Val_Cell_Year)<-TimeNames
    #get the sum of cell values for each year for the species
    #put them in the compilation data frame
    SpYearSumm[mySpp,]<-colSums(Spp_Val_Cell_Year,na.rm=T)
    print(toc())
    gc()
    if(writeSpRasters=="Yes"){
      for(myYear in as.character(writeYears))
        if (sp%in%writeSp|is.null(writeSp)){
          emptySpraster <- r
          values(r) <- Spp_Val_Cell_Year[,myYear]*100
          writeRaster(emptySpraster,
            OutTif,
            options = c("COMPRESS=LZW", "TFW=YES"),
            datatype = 'INT1U',
            overwrite = TRUE
          )
        }
      }
    
    
  }
    
    
    
   
    
    
    
 ######need to data   spyear sum ( currnetly in wide form) for below
  return(SpYearSumm)
}
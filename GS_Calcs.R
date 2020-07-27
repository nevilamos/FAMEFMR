##makeGS_LU------------
#Function makes LU matrix for Growth stage from TSF and EFG 
# YSF has 1 added to both the Lookup and the input to deal with YSF==0 which cannot be used in the array indexing
makeGS_LU<-function(EFG_TSF_4GS =myEFG_TSF_4GS){
  y<-EFG_TSF_4GS
  b=(y$YSF)+1
  c=y$EFG_NO
  e=y$GS4_NO
  x<-array(NA,dim=c(max(b),40))
  for(j in 1:nrow(y)){
    x[b[j],c[j]]<-e[j]
  }
  return(x)
}

 #############################################################
makeGS_Sum<-function(TimeSpan = FHanalysis$TimeSpan,
                     writeGSRasters,
                     myLU = GS_LU,
                     myResultsDir = ResultsDir,
                     myCropDetails = cropRasters,
                     myFHResults = FHanalysis,
                     myYSF_LFT = tsf_ysf_mat,
                     writeYears=NULL) {
  TimeRange<-as.integer(FHanalysis$TimeSpan)
  TimeNames<-as.character(FHanalysis$TimeSpan)
  LTR<-length(TimeRange)
  r<-FHanalysis$FH_IDr
  FH_ID<-values(r)
  r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
  gc()
  
  names(FH_ID)<-"ID"
  
  PLM<-cropRasters$PLM
  EFG<-cropRasters$EFG
  FIRE_REG<-as.integer(cropRasters$RGN)
  FIREFMZ<-as.integer(cropRasters$FIREFMZ)
  DELWP<-cropRasters$DELWP
  #PUBLIC<-as.integer()
  AllCombs<-as.data.table(cbind(FH_ID,EFG,FIRE_REG,FIREFMZ,PLM,DELWP))
  #TFI<-left_join(EFG_DF,TFI_LUT)
  rm(PLM,EFG,FIRE_REG,FIREFMZ,FH_ID,DELWP)
  gc()
  #calc the unique combinations of fire history(FH_ID) and EFG and various other admin datasets-----
  #This then allows the reduction of the subsequent calculation matrices to the number of unique cases of EFG and fire history rather than the numebr of raster cells
  U_AllCombs<-mgcv::uniquecombs(AllCombs)
  #index of the unique combs for linking back to raster
  Index_AllCombs<-attributes(U_AllCombs)$index
  rm(AllCombs)
  gc()
  
  
  #get the number of pixels in each this can be used at the end of the process to calcualte area cases
  nPixel<-as.vector(table(Index_AllCombs))
  #Add index column a and count of pixels to unique combs matrix
  Index<-1:nrow(U_AllCombs)
  U_AllCombs<-cbind(Index,nPixel,U_AllCombs)
  
  
  
  ### using DT formatted left join ( the "left" table is the one in brackets on the right)
  setDT(TFI_LUT)
  setkey(TFI_LUT,"EFG")
  U_AllCombs<-as.data.table(U_AllCombs)
  setkey(U_AllCombs,"EFG")
  U_AllCombs_TFI<-TFI_LUT[U_AllCombs]
  #have to reset the index to Index to retrun to the original order which is needed for cbinds below.
  setkey(U_AllCombs_TFI,"Index")
  
  
  #read the corresponding shapefile for the FireHat analysis chosen
  #convert to a data.frame
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  ID<-OutTab$ID
  
  YSF_Fields<-names(OutTab)[grep("^YSF",names(OutTab))]
  YSFplus1<-OutTab[,YSF_Fields]+1
  #partial inflation using U_AllCombs_TFI$FH_ID
  YSFplus1<-YSFplus1[U_AllCombs_TFI$FH_ID,]
  YSFplus1<-as.matrix(YSFplus1)
  #inflate GS_LU on EFG axis using U_AllCombs_TFI$$EFG
  GS_LU<-t(GS_LU)[U_AllCombs_TFI$EFG,]
  
  #need to work out how to use the TSF+1 value to identify the column value 
  Nrow<-nrow(YSFplus1)
  Ncol<-ncol(YSFplus1)
  GSrows<-matrix(,Nrow,Ncol)
  for(i in 1:Nrow){
    GSrows[i,]<-GS_LU[i,YSFplus1[i,]]
  }
  
 ###################################  
 ################################### 
  
  colnames(GSYrArray)<-TimeSpan
  myTab<-cbind(EFG,FIRE_REG,DELWP,FIREFMZ,PLM,GSYrArray)
  myTab<-as_tibble(myTab)
  #tablong<-gather(myTab,SEASON,GS,-EFG,-FIRE_REG,-FIREFMZ,-PLM)
  GS_Summary<-myTab%>%
    gather(SEASON,GS,-EFG,-FIRE_REG,-FIREFMZ,-PLM,-DELWP)%>%
    count(EFG,FIRE_REG,DELWP,FIREFMZ,PLM,SEASON,GS)
  rm(myTab)
  gc()
  GS_Summary$Hectares<-GS_Summary$n*cellsToHectares()
  GS_Summary<-left_join(GS_Summary,GS_LUT)
  GS_Summary<-left_join(GS_Summary,TFI_LUT[,c("EFG","EFG_NAME")])
  GS_Summary<-left_join(GS_Summary,FIREFMZ_LUT)
  GS_Summary<-left_join(GS_Summary,REG_LUT)
  GS_Summary<-left_join(GS_Summary,DELWP_LUT)
  return(GS_Summary)
}

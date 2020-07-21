
LBY_f<-function(M,y){
  M[M>y|M==0]<-NA
  LBY<-Rfast::rowMaxs(M,value = T)
  LBY[is.infinite(LBY)]<-NA
  return(LBY)
}
#Calcuates where each cell is currently at below MinTFI or above HI_TFI------------------------------------
#returns the per cell and long and wide tablessummariesed by season and EFG
calc_TFI_2<-function(FHanalysis,
                     cropRasters,
                     TFI_LUT,
                     OutputRasters = makeTFIRasters){
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
  #PUBLIC<-as.integer()
  AllCombs<-as.data.table(cbind(FH_ID,EFG,FIRE_REG,FIREFMZ,PLM))
  #TFI<-left_join(EFG_DF,TFI_LUT)
  rm(PLM,EFG,FIRE_REG,FIREFMZ,FH_ID)
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
  #INTFields<-names(OutTab)[grep("^INT",names(OutTab))]
  SEASFields<-names(OutTab)[grep("^SEAS",names(OutTab))]
  TYPEFields<-names(OutTab)[grep("^FireType",names(OutTab))]
  
  
  SEAS<-as.matrix(OutTab[,SEASFields])#The Season
  SEAS[SEAS==0]<-NA
  TYPE<-as.matrix(OutTab[,TYPEFields[]])#The type of the fire 
  
  TYPE_HI<-TYPE==2
  TYPE_HI[TYPE_HI<1]<-NA
  SEAS_HI<-SEAS*TYPE_HI
  
  TYPE_LO<-TYPE==1
  TYPE_LO[TYPE_LO<1]<-NA
  SEAS_LO<-SEAS*TYPE_LO
  rm(OutTab)
  gc()
  
  
  
  # Calc Last Burned Year (LBY) for each firetype (currently only deals with two firetypes)
  LBY_HI<-matrix(NA,nrow(SEAS),LTR)
  colnames(LBY_HI)<-TimeNames
  LBY_HI<-LBY_LO<-cbind(LBY_HI,ID)
  
  for(i in 1:LTR){
    try({
      y=TimeRange[i]
      
      LBY_HI[,i]<-LBY_f(M=SEAS_HI,y)
      LBY_LO[,i]<-LBY_f(M=SEAS,y)##############This Should Maybe be SEAS_LO
      print(y)
      
    })
  }
  
  
  
  
  #partial inflation using U_AllCombs_TFI$FH_ID
  #and TFI vaules using U_AllCombs_TFI$MIN_LO_TFI
  LBY_LO<-LBY_LO[U_AllCombs_TFI$FH_ID,]
  LBY_HI<-LBY_HI[U_AllCombs_TFI$FH_ID,]
  
  
  ###Calc the TFI status - this is the section to check if there are unusual TFI statuses -------
  TFI_LO<-(t(TimeRange-t(LBY_LO[,TimeNames]))-U_AllCombs_TFI$MIN_LO_TFI)<1
  TFI_HI<-(t(TimeRange-t(LBY_HI[,TimeNames]))-U_AllCombs_TFI$MIN_HI_TFI)<1
  TFI_MAX<-(t(TimeRange-t(LBY_LO[,TimeNames]))>U_AllCombs_TFI$MAX_TFI)*5L
  TFI_COMB<-TFI_LO
  TFI_COMB[TFI_HI==T]<-TRUE
  #colnames(TFI_COMB)
  TFI_STATUS<-TFI_COMB+TFI_MAX
  ###Turn next line on to get rid of TFI status>2----
  #TFI_STATUS[TFI_STATUS>4]<-2
  TFI_STATUS[is.na(TFI_STATUS)]<--99L
  #TFI_STATUS<-cbind(TFI_STATUS,U_AllCombs_TFI)
  rm()
  
  #dplyr wrangling of output summary tables-------------  
  TFI_Long<-cbind(TFI_STATUS,U_AllCombs_TFI[,c("EFG",  "EFG_NAME", "nPixel",  "FIRE_REG", "FIREFMZ", "PLM")])%>%
    pivot_longer(all_of(TimeNames),names_to="SEASON",values_to="TFI_STATUS")%>%
    group_by(EFG,  EFG_NAME,  FIRE_REG, FIREFMZ, PLM,SEASON,TFI_STATUS)%>%
    summarize(Cells=sum(nPixel))%>%
    mutate(hectares=cellsToHectares(Cells))
  
  # %>%
  #   mutate(hectares=(nPixel))%>%
  #   select(-nPixel)
  # 
  # 
  # TFI_Summary<-as_tibble(z%>%
  #                          gather(SEASON,TFI_VAL,-c(EFG,FIRE_REG,FIREFMZ,PLM))%>%
  #                          count(EFG,FIRE_REG,FIREFMZ,PLM,SEASON,TFI_VAL))
  # 
  # 
  # TFI_Summary<-left_join(TFI_Summary,TFI_STATUS_LUT)
  # TFI_Summary<-left_join(TFI_Summary,TFI_LUT[,c("EFG","EFG_NAME")])
  # TFI_Summary<-left_join(FIREFMZ_LUT,TFI_Summary)
  # TFI_Summary<-left_join(TFI_Summary,REG_LUT)
  
  #expand the TFI status values to raster vector for raster values save as filematrix - for future reference or reading into raster. --------------
  TFI_STATUS_RASTER_VALS<-fm.create(file.path(ResultsDir,"TFI_STATUS_RASTER_VALS"),nrow=length(Index_AllCombs),ncol=LTR)
  colnames(TFI_STATUS_RASTER_VALS)<-TimeNames
  for(i in 1:LTR){
    j=TimeNames[i]
    TFI_STATUS_RASTER_VALS[,i]<-TFI_STATUS[Index_AllCombs,j]
    print(j)
  }
  #if filematrix of raster values already made in order to open again----------------
  TFI_STATUS_RASTER_VALS<-fm.open(file.path(ResultsDir,"TFI_STATUS_RASTER_VALS"))
  
  # write the TFI Status tiffs for each year.-------------- 
  #should look at finding way to assign to stack or even to envi or similar directly - maybe talk to Pete Griff and Lachlan about this.
  
  if (OutputRasters == "Yes"){
    #cl<-makeCluster(Ncores,outfile="")
    #registerDoParallel(cl, cores=Ncores)
    #foreach(i=iter(1:LTR),.packages =c("raster","filematrix") )%dopar%{
    
    for (i in 1:LTR){
      j=TimeNames[i]
      outR<-r
      values(outR)<-TFI_STATUS_RASTER_VALS[,i]
      writeRaster(outR,file.path(ResultsDir,"TFI_Rasters",paste0("TFI_STATUS",j,".tif")),overwrite=T,datatype="INT1U")
      rm(outR)
      gc()
    }
    #stopImplicitCluster()
    
  }
  return(TFI_Long)
}
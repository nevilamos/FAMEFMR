calc_TFI_75<-function(FHanalysis=FHanalysis,
                   cropRasters=cropRasters,
                   TFI_LUT_DF = TFI_LUT,
                   OutputRasters = makeTFIRasters)#whether or not to output rasters for each year of area under TFI
{
  TimeRange<-as.integer(FHanalysis$TimeSpan)
  if (Test){TimeRange=TimeRange[1:2]}
  r<-FHanalysis$FH_IDr
#FH_ID<-as.data.frame(as.integer(values(r)))
#names(FH_ID)<-"ID"
  FH_ID<-as.integer(values(r))
  
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  OutTab<-as.data.table(OutTab)
  rm(FHanalysis)
  gc()
  
  #read in the raster of EFG numbers PLM FireFMZ etc
  #extract the values as a vector and use it to make a data frame contianing cell-wise values for TFI
  # PLM<-cropRasters$PLM
  # EFG<-cropRasters$EFG
  # FIRE_REG<-as.integer(cropRasters$RGN)
  # FIREFMZ<-as.integer(cropRasters$FIREFMZ)
  #PUBLIC<-as.integer()
  EFG_DF<-data.table(cropRasters$EFG,cropRasters$RGN,cropRasters$FIREFMZ,cropRasters$PLM)#cbind(EFG,FIRE_REG,FIREFMZ,PLM)
  rm(cropRasters)
  gc()
  names(EFG_DF)<-c("EFG","FIRE_REG","FIREFMZ","PLM")
  TFI<-as.data.table(TFI_LUT_DF[,-5])[EFG_DF, on="EFG"]   
  #TFI<-left_join(EFG_DF,TFI_LUT_DF[,-4])
  gc()
  #read the corresponding shapefile for the FireHat analysis chosen
  #convert to a data.frame
  
 
    SEASFields<-names(OutTab)[grep("^SEAS",names(OutTab))]
  TYPEFields<-names(OutTab)[grep("^FireType",names(OutTab))]
  outM<-as.matrix(OutTab)
  mode(outM)<-"integer"
  SEAS<-outM[FH_ID,SEASFields]
  #INTFields<-names(OutTab)[grep("^INT",names(OutTab))]

  gc()
  #using data.table notation for these ( left_joins) as faster
  #SEAS1<-as.matrix(OutTab[FH_ID,..SEASFields])
  mode(SEAS)<-"integer"
  #SEAS<-as.matrix(x[,SEASFields])#The Season
  SEAS[SEAS==0L]<-NA
  gc()
  TYPE<-outM[FH_ID,TYPEFields]
  mode(TYPE)<-"integer"
  gc()
  #SEAS_HI<-x<-left_join(FH_ID,OutTab)
  #x<-left_join(FH_ID,OutTab)

  
  #SEAS<-as.matrix(x[,SEASFields])#The Season

  #TYPE<-as.matrix(x[,TYPEFields[]])#The type of the fire 
  TYPE_HI<-TYPE==2L
  TYPE_HI[TYPE_HI<1]<-NA
  gc()
  TYPE_LO<-TYPE==1L
  gc()
  SEAS_HI<-SEAS*TYPE_HI
  
  #rm(x)
  gc()
  #print("making multicore cluster for parallel processing of TSF")
  #Ncores<-8#detectCores()-2
  #cl<-makeCluster(Ncores,outfile="")
  #registerDoParallel(cl, cores=Ncores)
  
  LBY_HI<-NULL
  #LBY_HI<-foreach(y=iter(TimeRange),.combine = cbind,.packages ="Rfast" )%dopar%
  for(y in TimeRange){
    try({
      X<-SEAS_HI
      X[X>y]<-NA
      gc()
      X[X==0]<-NA
      mode(X)<-"numeric"
      #LBY<-apply(X,1,max,na.rm=TRUE)
      LBY<-rowMaxs(X,value=T)
      rm(X)
      gc()
      LBY[is.infinite(LBY)]<-NA
      mode(LBY)<-"integer"
      LBY
      LBY_HI<-cbind(LBY_HI,LBY)
      print(y)
      gc()
    })
  }
  mode(LBY_HI)<-"integer"
  gc()
  #LBY_LO<-foreach(y=iter(TimeRange),.combine = cbind,.packages ="Rfast" )%dopar%
  LBY_LO<-NULL  
  for(y in TimeRange){
    try({
      X<-SEAS
      X[X>y]<-NA
      X[X==0]<-NA
      mode(X)<-"numeric"
      LBY<-rowMaxs(X,value=T)
      rm(X)
      gc()
      LBY[is.infinite(LBY)]<-NA
      LBY_LO<-cbind(LBY_LO,LBY)
      gc()
    })
  }
  mode(LBY_LO)<-"integer"
  
  TFI_LO<-(t(TimeRange-t(LBY_LO))-TFI$MIN_LO_TFI)<1
  TFI_HI<-(t(TimeRange-t(LBY_HI))-TFI$MIN_HI_TFI)<1
  TFI_MAX<-(t(TimeRange-t(LBY_LO))>TFI$MAX_TFI)*2L
  TFI_COMB<-TFI_LO
  TFI_COMB[TFI_HI==T]<-TRUE
  colnames(TFI_COMB)<-as.character(TimeRange)
  TFI_STATUS<-TFI_COMB+TFI_MAX
  TFI_STATUS[is.na(TFI_STATUS)]<--99L
  #Area_ha<-(as.numeric(FHanalysis$RasterRes)/100)^2
  z<-as_tibble(cbind(EFG_DF,TFI_STATUS))
  
  #zz<-gather(z[,!names(z)%in%c("MIN_LO_TFI", "MIN_HI_TFI", "MAX_TFI", "EFG_NAME")]
  TFI_Summary<-as_tibble(z%>%
                           gather(SEASON,TFI_VAL,-c(EFG,FIRE_REG,FIREFMZ,PLM))%>%
                           count(EFG,FIRE_REG,FIREFMZ,PLM,SEASON,TFI_VAL))
  TFI_Summary$Hectares<-TFI_Summary$n*cellsToHectares()
  
  TFI_Summary<-left_join(TFI_Summary,TFI_STATUS_LUT)
  TFI_Summary<-left_join(TFI_Summary,TFI_LUT[,c("EFG","EFG_NAME")])
  TFI_Summary<-left_join(FIREFMZ_LUT,TFI_Summary)
  TFI_Summary<-left_join(TFI_Summary,REG_LUT)
  
  
  
  if (OutputRasters == "Yes"){
    #write a lookup table for the values in the TFI _STATUS_TIFs
    
    # write the TFI Status tiffs for each year.  
    outR<-r
    for(i in colnames(TFI_STATUS)){
      values(r)<-TFI_STATUS[,i]
      writeRaster(r,file.path(ResultsDir,"TFI_Rasters",paste0("TFI_STATUS",i,".tif")),overwrite=T,datatype="INT1U")
    gc()
    }
    
  }
  
  return(TFI_Summary)
  
  }
  
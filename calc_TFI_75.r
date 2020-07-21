calc_TFI_75<-function(FHanalysis=FHanalysis,
                      cropRasters=cropRasters,
                      TFI_LUT_DF = TFI_LUT,
                      OutputRasters = makeTFIRasters)#whether or not to output rasters for each year of area under TFI
{
  TimeRange<-as.integer(FHanalysis$TimeSpan)
  #if (Test){TimeRange=TimeRange[1:2]}
  r<-FHanalysis$FH_IDr
  Ncells<-ncell(r)
  NSeas<-length(TimeRange)
  
  FH_ID<-as.integer(values(r))
  
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  OutTab<-as.data.table(OutTab)
  #rm(FHanalysis)
  gc()
  
  #read in the raster of EFG numbers PLM FireFMZ etc
  #extract the values as a vector and use it to make a data frame contianing cell-wise values for TFI
  
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
  TFI_STATUS<-fm.create(filenamebase = "TFI_STATUS",nrow=Ncells,ncol=NSeas,type="integer",)
  for(i in 1:NSeas){
    try({
      y=TimeRange[i]
      
      X<-SEAS_HI
      X[X>y|X==0]<-NA
      #X[X==0]<-NA
      mode(X)<-"numeric"
      
      LBY_HI<-rowMaxs(X,value=T)
      rm(X)
      gc()
      LBY_HI[is.infinite(LBY_HI)]<-NA
      mode(LBY_HI)<-"integer"
      
      X<-SEAS
      X[X>y]<-NA
      X[X==0]<-NA
      mode(X)<-"numeric"
      LBY_LO<-rowMaxs(X,value=T)
      rm(X)
      gc()
      LBY_LO[is.infinite(LBY_LO)]<-NA
      
      TFI_LO<-(TimeRange-LBY_LO-TFI$MIN_LO_TFI)<1
      TFI_HI<-(TimeRange-LBY_HI-TFI$MIN_HI_TFI)<1
      TFI_MAX<-((TimeRange-LBY_LO)>TFI$MAX_TFI)*2L
      TFI_COMB<-TFI_LO
      TFI_COMB[TFI_HI==T]<-TRUE
      TFI_STAT<-TFI_COMB+TFI_MAX
      TFI_STAT[is.na(TFI_STAT)]<--99L
      TFI_STATUS[,i]<-TFI_STAT
      rm()
      print(y)
    })
  }
  TFI_STATUS_M<-as.matrix(TFI_STATUS)
  colnames(TFI_STATUS_M)<-as.character(TimeRange)
  
  
  z<-as_tibble(cbind(EFG_DF,TFI_STATUS_M))
  
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
    for(i in 1:NSeas){
      values(r)<-TFI_STATUS[,i]
      writeRaster(r,file.path(ResultsDir,"TFI_Rasters",paste0("TFI_STATUS",TimeRange[i],".tif")),overwrite=T,datatype="INT1U")
      gc()
    }
    
  }
  
  return(TFI_Summary)
  
}

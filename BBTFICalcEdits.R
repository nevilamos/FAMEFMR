calc_BBTFI<-function(FHanalysis,#the slected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                     cropRasters,
                     TFI_LUT_DF = TFI_LUT#the dataframe read from a csv that gives the lookup table from EVD to MinTFI_LO, MinTFI_HI and MaxTFI
) {
  
  r<-FHanalysis$FH_IDr
  Hectares<-(as.numeric(FHanalysis$RasterRes)/100)^2
  
  FH_ID<-as.data.frame(values(r))
  names(FH_ID)<-"ID"
  PLM<-cropRasters$PLM
  EFG<-cropRasters$EFG
  FIRE_REG<-cropRasters$RGN
  FIREFMZ<-cropRasters$FIREFMZ
  PLM<-cropRasters$PLM
  EFG_DF<-as.data.frame(EFG)
  
  TFI<-left_join(EFG_DF,TFI_LUT_DF)
  
  #read the corresponding shapefile for the FireHat analysis chosen
  #convert to a data.frame
  
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  
  INTFields<-names(OutTab)[grep("^INT",names(OutTab))]#fields from FH analysis giving inter fire intervalsof 1:nth fire
  SEASFields<-names(OutTab)[grep("^SEAS",names(OutTab))]#fields giving season of 1:nth fire
  TYPEFields<-names(OutTab)[grep("^FireType",names(OutTab))]#fields giving type of 1:nth fire
  x<-left_join(FH_ID,OutTab)
  x<-as.matrix(x)
  mode(x)<-"integer"
  INTS<-x[,INTFields]#The fire intervals in numbered sequence
  SEAS<-x[,SEASFields[-1]]#The Season of the second fire of each interval to get the season when burning under TFI occurs
  SEAS[SEAS==0]<-NA
  TYPE<-x[,TYPEFields[-1*length(TYPEFields)]]#The type of the first fire of the interval to determine whether high or lo TFI applies
  
  
  TYPE2<-as.matrix(x[,TYPEFields[-1]])#The firetype of the second fire - that determines whether the fire at the date that the second burn event occured was high or low for reporting
  TYPE2[TYPE2==0]<-NA
  
  
  rm(x)
  gc()
  
  TYPE_HI<-TYPE==2
  TYPE_LO<-TYPE==1
  
  
  
  
  BB_LO_TFI_SEASON<-SEAS*(INTS<TFI$MIN_LO_TFI)#if any fire interval is less than the MIN_LO_TFI then the veg has been burnt below TFI
  BB_LO_TFI_SEASON[BB_LO_TFI_SEASON==0]<-NA#remove zeros because later going to want to calc minimum dates
  BB_HI_TFI_SEASON<-SEAS*TYPE_HI*(INTS<TFI$MIN_HI_TFI)#only in cases where th first firetype is HI and the interval is below MIN_HI_TFI is the veg burnt below TFI
  BB_HI_TFI_SEASON[BB_HI_TFI_SEASON==0]<-NA
  
  #next two line combine dates for fires BBTFI that are below the low TFI threhsold,
  #and those below the Hi TFI theshold( as determined by the TYPE of the first fire)
  BBTFI_COMB<-BB_LO_TFI_SEASON
  BBTFI_COMB[is.na(BB_LO_TFI_SEASON)]<-BB_HI_TFI_SEASON[is.na(BB_LO_TFI_SEASON)]
  
  #separating whether the fire in the second season of the BBTFI was high or low for reporting
  ID<-1:nrow(BBTFI_COMB)
  BBTFI_ID<-as_tibble(cbind(ID,BBTFI_COMB))%>%gather(Garbage,SEASON,-ID)
  TYPE2_ID<-as_tibble(cbind(ID,TYPE2))%>%gather(Garbage,TYPE,-ID)
  
  ID_LU<-as_tibble(cbind(ID,EFG,FIREFMZ,FIRE_REG,PLM))
  
  BBTFI_BY_TYPE<-right_join(ID_LU,na.omit(cbind(BBTFI_ID[,c("ID","SEASON")],TYPE2_ID[,"TYPE"])))%>%
    count(EFG,FIRE_REG,FIREFMZ,PLM,SEASON,TYPE)
  BBTFI_BY_TYPE$Hectares<-BBTFI_BY_TYPE$n*Hectares
  BBTFI_BY_TYPE<-right_join(FIRETYPE_LUT,BBTFI_BY_TYPE)
  BBTFI_BY_TYPE<-right_join(TFI_LUT[,c("EFG","EFG_NAME")],BBTFI_BY_TYPE)
  BBTFI_BY_TYPE<-right_join(FIREFMZ_LUT,BBTFI_BY_TYPE)
  BBTFI_BY_TYPE<-right_join(REG_LUT,BBTFI_BY_TYPE)
  
  x<-BBTFI_COMB
  EFGM<-matrix(EFG,dim(x)[1],dim(x)[2])
  y<-!is.na(x)
  z<-t(apply(y,1,cumsum))
  yy<-z[y]
  zz<-EFGM[y]
  xx<-x[y]
  BBTFI_Cell_SEASON<-as.data.frame(cbind(zz,xx,yy,Hectares))
  names(BBTFI_Cell_SEASON)[1:3]<-c("EFG","SEASON","Times_BBTFI")
  BBTFI_EFG_Area_SEASON<-BBTFI_Cell_SEASON%>%group_by(EFG,SEASON,Times_BBTFI)%>%summarise(ha=sum(Hectares))
  #calculate the number of times each  row ( cell of matrix) has been BBTFI using apply sum on true/false
  #TimesBBTFI<-apply(!is.na(BBTFI_COMB),1, sum)
  TimesBBTFI<-rowSums(!is.na(BBTFI_COMB))
  FirstBBTFI<-apply(BBTFI_COMB,1,min,na.rm=T )
  FirstBBTFI[is.infinite(FirstBBTFI)]<-NA
  BBTFIRaster<-cropRasters$Raster
  values(BBTFIRaster)<-TimesBBTFI
  writeRaster(BBTFIRaster,file.path(ResultsDir,"TFI_Rasters","TimesBBTFI.tif"),overwrite=T)
  values(BBTFIRaster)<-FirstBBTFI
  writeRaster(BBTFIRaster,file.path(ResultsDir,"TFI_Rasters","FirstBBTFI.tif"),overwrite=T)
  y<-data.frame(EFG,TimesBBTFI)
  y$Area_ha<-Hectares
  yy<-y%>%group_by(EFG,TimesBBTFI)%>%summarise(ha=sum(Area_ha))
  TimesBBTFI_Summary<-dcast(yy,formula =EFG~TimesBBTFI,value.var = "ha")
  
  #apply(myBBTFI$BBTFI_COMB,1,min,na.rm=T)
  
  
  #BBTFI_COMB<-list("BBTFI_COMB"=BBTFI_COMB,"TimesBBTFI_Summary"=TimesBBTFI_Summary,"BBTFI_Cell_SEASON"=BBTFI_Cell_SEASON,"BBTFI_EFG_Area_SEASON"=BBTFI_EFG_Area_SEASON,"TFI"=TFI,"BBTFI_BY_TYPE"=BBTFI_BY_TYPE)
  write.csv(TimesBBTFI_Summary,file.path(ResultsDir,"TimesBBTFI_Summary.csv"))
  #write.csv(BBTFI_COMB,file.path(ResultsDir,"BBTFI_COMB.csv"))
  write.csv(BBTFI_EFG_Area_SEASON,file.path(ResultsDir,"BBTFI_EFG_Area_SEASON.csv"))
  BBTFIresults<-list("BBTFI_COMB"=BBTFI_COMB,"TimesBBTFI_Summary"=TimesBBTFI_Summary,"BBTFI_Cell_SEASON"=BBTFI_Cell_SEASON,"BBTFI_EFG_Area_SEASON"=BBTFI_EFG_Area_SEASON,"TFI"=TFI,"BBTFI_BY_TYPE"=BBTFI_BY_TYPE)
  return(BBTFIresults)
}

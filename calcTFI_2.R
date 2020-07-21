##Calcuates where each cell is currently at below MinTFI or above HI_TFI------------------------------------
#returns the per cell and long and wide tablessummariesed by season and EFG
calc_TFI2<-function(FHanalysis,
                    cropRasters,
                    TFI_LUT_DF = TFI_LUT,
                    OutputRasters = makeTFIRasters)#whether or not to output rasters for each year of area under TFI
{
  TimeRange<-as.integer(FHanalysis$TimeSpan)
  r<-FHanalysis$FH_IDr
  FH_ID<-values(r)
  #FH_ID<-as.data.frame(as.integer(values(r)))
  names(FH_ID)<-"ID"
  
  #read in the raster of EFG numbers PLM FireFMZ etc
  #extract the values as a vector and use it to make a data frame contianing cell-wise values for TFI
  PLM<-cropRasters$PLM
  EFG<-cropRasters$EFG
  FIRE_REG<-as.integer(cropRasters$RGN)
  FIREFMZ<-as.integer(cropRasters$FIREFMZ)
  #PUBLIC<-as.integer()
  EFG_DF<-as.data.frame(cbind(EFG,FIRE_REG,FIREFMZ,PLM))
  #TFI<-left_join(EFG_DF,TFI_LUT_DF)
  ID_EFG<-cbind(FH_ID,EFG)
  rm(PLM,EFG,FIRE_REG,FIREFMZ)
  gc()
  #calc the unique combinations of fire history(FH_ID) and EFG
  U_ID_EFG<-mgcv::uniquecombs(ID_EFG)
  Index<-1:nrow(U_ID_EFG)
  U_ID_EFG<-cbind(Index,U_ID_EFG)
  Index_U_ID_EFG<-attributes(U_ID_EFG)$index
  
  
  setDT(TFI_LUT_DF)
  setkey(TFI_LUT_DF,"EFG")
  U_ID_EFG<-as.data.table(U_ID_EFG)
  setkey(U_ID_EFG,"EFG")
  U_ID_EFG_TFI<-TFI_LUT_DF[U_ID_EFG]

  
  #read the corresponding shapefile for the FireHat analysis chosen
  #convert to a data.frame
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  ID<-OutTab$ID
  #INTFields<-names(OutTab)[grep("^INT",names(OutTab))]
  SEASFields<-names(OutTab)[grep("^SEAS",names(OutTab))]
  TYPEFields<-names(OutTab)[grep("^FireType",names(OutTab))]
  #SEAS_HI<-x<-left_join(FH_ID,OutTab)
  
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
  
  #print("making multicore cluster for parallel processing of TSF")
  #Ncores<-8#detectCores()-2
  #cl<-makeCluster(Ncores,outfile="")
  #registerDoParallel(cl, cores=Ncores)
  
  LBY_HI<-matrix(NA,nrow(SEAS),length(TimeRange))
  colnames(LBY_HI)<-as.character(TimeRange)
  LBY_HI<-LBY_LO<-cbind(LBY_HI,ID)
  
  #LBY_HI<-foreach(y=iter(TimeRange),.combine = cbind,.packages ="Rfast" )%dopar%
  for(i in 1:length(TimeRange)){
    try({
      y=TimeRange[i]
      
      LBY_HI[,i]<-LBY_f(M=SEAS_HI,y)
      
      print(y)
      
    })
  }
  
  for(i in 1:length(TimeRange)){
    try({
      y=TimeRange[i]
      
      LBY_LO[,i]<-LBY_f(M=SEAS,y)##############This ShOULD MayBE be SEAS_LO
      
      print(y)
      
    })
  }
 
  ##########the stuff above needs to be inflated at some point below to match the raster cells 
  
  #infalte using U_ID_EFG_TFI$FH_ID
  #and TFI vales using U_ID_EFG_TFI$MIN_LO_TFI
  
  
  
  
  
  TFI_LO<-(t(TimeRange-t(LBY_LO))-TFI$MIN_LO_TFI)<1
  TFI_HI<-(t(TimeRange-t(LBY_HI))-TFI$MIN_HI_TFI)<1
  TFI_MAX<-(t(TimeRange-t(LBY_LO))>TFI$MAX_TFI)*2L
  TFI_COMB<-TFI_LO
  TFI_COMB[TFI_HI==T]<-TRUE
  colnames(TFI_COMB)<-as.character(TimeRange)
  TFI_STATUS<-TFI_COMB+TFI_MAX
  TFI_STATUS[is.na(TFI_STATUS)]<--99L
  Area_ha<-(as.numeric(FHanalysis$RasterRes)/100)^2
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
    }
    
  }
  
  return(TFI_Summary)
  
}
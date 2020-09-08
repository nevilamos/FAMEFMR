# Function to join ominated Lookup Tables to data frame containing indces of LUT
Join_Names<-function(myDF,#dataframe or similar containing indices for the LUTS listed
                     LUTS=c("TFI_LUT","TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")){
  for(i in LUTS){
    try(myDF<-left_join(myDF,get(i)))}
  return (myDF)
}


# Function to calculate last burnt year from matrix of rows of fire season  iterating by year (y) used in calc_TFI_2---------
LBY_f<-function(M,y){
  M[M>y|M==0]<-NA
  LBY<-Rfast::rowMaxs(M,value = T)
  LBY[is.infinite(LBY)]<-NA
  return(LBY)
}



#Calcuates where each cell is currently at below MinTFI or above MAX_TFI------------------------------------
#returns the per cell and long table summarised by multiple admin units and evc
calc_TFI_2<-function(FHanalysis,
                     U_AllCombs_TFI=myAllCombs$U_AllCombs_TFI,
                     Index_AllCombs=myAllCombs$Index_AllCombs,
                     TFI_LUT,
                     OutputRasters = makeTFIRasters){
  
  
  TimeRange<-as.integer(FHanalysis$TimeSpan)
  TimeNames<-as.character(FHanalysis$TimeSpan)
  LTR<-length(TimeRange)
  r<-FHanalysis$FH_IDr
  FH_ID<-values(r)
  r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
  gc()

  
  
  #read the corresponding shapefile for the FireHat analysis chosen
  #convert to a data.frame
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  ID<-OutTab$ID
  #INTFields<-names(OutTab)[grep("^INT",names(OutTab))]
  SEASFields<-names(OutTab)[grep("^SEAS",names(OutTab))]
  TYPEFields<-names(OutTab)[grep("^FireType",names(OutTab))]
  
  
  SEAS<-as.matrix(OutTab[,SEASFields])#The Season
  SEAS[is.na(SEAS)]<-0
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
  TFI_VAL<-TFI_COMB+TFI_MAX
  ###Turn next line on to get rid of TFI status>2----
  #TFI_VAL[TFI_VAL>4]<-2
  TFI_VAL[is.na(TFI_VAL)]<--99L
  #TFI_VAL<-cbind(TFI_VAL,U_AllCombs_TFI)
  
  #This next section is only run for debgging unusual TFI statuses it allows their isolation at the level of unique combination of fire history and EFG ----
  # Check_TFI<-cbind(TFI_VAL,U_AllCombs_TFI)%>%
  #   select(-FIRE_REG,-FIREFMZ,-PLM,-DELWP)%>%
  #   pivot_longer(all_of(TimeNames),names_to="SEASON",values_to="TFI_VAL")%>%
  #   filter(!TFI_VAL%in%c(-99,0,1,5))%>%
  #   group_by(EFG,MIN_LO_TFI,MIN_HI_TFI,MAX_TFI,EFG_NAME,Index,FH_ID, SEASON,TFI_VAL)%>%
  #   summarize(Cells=sum(nPixel))
  # 
  # Check_TFI<-as.data.table(Check_TFI)
  # setkey(Check_TFI,"FH_ID")
  # write.csv(Check_TFI,"Check_TFI.csv")
  # OutTab<-as.data.table(OutTab)
  # Check_TFI<-OutTab[Check_TFI]
  
  #dplyr wrangling of output summary tables-------------  
  TFI_Summary<-cbind(TFI_VAL,U_AllCombs_TFI)%>%
    pivot_longer(all_of(TimeNames),names_to="SEASON",values_to="TFI_VAL")%>%
    group_by(EFG_NAME,FIRE_FMZ_NAME,FIRE_FMZ_SHORT_NAME,	FIRE_REGION_NAME,	DELWP_REGION,	EFG,	FIRE_REG,	FIREFMZ,	PLM,	DELWP,SEASON,TFI_VAL)%>%
    summarize(nCells=sum(nPixel),Hectares=sum(Hectares))
  
  TFI_Summary<-left_join(TFI_Summary,TFI_STATUS_LUT)
  

  if (OutputRasters == "Yes"){
    # #expand the TFI status values to raster vector for raster values save as filematrix - for future reference or reading into raster. --------------
    # TFI_VAL_RASTER_VALS<-fm.create(file.path(ResultsDir,"TFI_VAL_RASTER_VALS"),nrow=length(Index_AllCombs),ncol=LTR)
    # colnames(TFI_VAL_RASTER_VALS)<-TimeNames
    # for(i in 1:LTR){
    #   j=TimeNames[i]
    #   TFI_VAL_RASTER_VALS[,i]<-TFI_VAL[Index_AllCombs,j]
    #   print(j)
    # }
    # #if filematrix of raster values already made in order to open again----------------
    #TFI_VAL_RASTER_VALS<-fm.open(file.path(ResultsDir,"TFI_VAL_RASTER_VALS"))
    
    # write the TFI Status tiffs for each year.-------------- 
    #should look at finding way to assign to stack or even to envi or similar directly - maybe talk to Pete Griff and Lachlan about this.
    
    
    #cl<-makeCluster(Ncores,outfile="")
    #registerDoParallel(cl, cores=Ncores)
    #foreach(i=iter(1:LTR),.packages =c("raster","filematrix") )%dopar%{
    
    for (i in 1:LTR){
      j=TimeNames[i]
      outR<-r
      values(outR)<-TFI_VAL[Index_AllCombs,j]
      writeRaster(outR,file.path(ResultsDir,"TFI_Rasters",paste0("TFI_VAL",j,".tif")),overwrite=T,datatype="INT1U")
      rm(outR)
      print(paste("Output TFI Raster", j))
      gc()
    }
    #stopImplicitCluster()
    
  }
  return(TFI_Summary)
}

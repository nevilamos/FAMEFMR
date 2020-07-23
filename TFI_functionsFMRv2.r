# Function to calculate last burnt year from matrix of rows of fire season  iterating by year (y) used in calc_TFI_2---------
LBY_f<-function(M,y){
  M[M>y|M==0]<-NA
  LBY<-Rfast::rowMaxs(M,value = T)
  LBY[is.infinite(LBY)]<-NA
  return(LBY)
}



#Calcuates where each cell is currently at below MinTFI or above HI_TFI------------------------------------
#returns the per cell and long table summarised by multiple admin units and evc
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
  TFI_Summary<-cbind(TFI_VAL,U_AllCombs_TFI[,c("EFG",  "EFG_NAME", "nPixel",  "FIRE_REG", "FIREFMZ", "PLM","DELWP")])%>%
    pivot_longer(all_of(TimeNames),names_to="SEASON",values_to="TFI_VAL")%>%
    group_by(EFG,  EFG_NAME,  FIRE_REG,DELWP, FIREFMZ, PLM,SEASON,TFI_VAL)%>%
    summarize(Cells=sum(nPixel))%>%
    mutate(hectares=cellsToHectares(Cells))
  
  TFI_Summary<-left_join(TFI_Summary,TFI_STATUS_LUT)
  
  TFI_Summary<-left_join(FIREFMZ_LUT,TFI_Summary)
  TFI_Summary<-left_join(TFI_Summary,REG_LUT)
  
  #expand the TFI status values to raster vector for raster values save as filematrix - for future reference or reading into raster. --------------
  TFI_VAL_RASTER_VALS<-fm.create(file.path(ResultsDir,"TFI_VAL_RASTER_VALS"),nrow=length(Index_AllCombs),ncol=LTR)
  colnames(TFI_VAL_RASTER_VALS)<-TimeNames
  for(i in 1:LTR){
    j=TimeNames[i]
    TFI_VAL_RASTER_VALS[,i]<-TFI_VAL[Index_AllCombs,j]
    print(j)
  }
  #if filematrix of raster values already made in order to open again----------------
  TFI_VAL_RASTER_VALS<-fm.open(file.path(ResultsDir,"TFI_VAL_RASTER_VALS"))
  
  # write the TFI Status tiffs for each year.-------------- 
  #should look at finding way to assign to stack or even to envi or similar directly - maybe talk to Pete Griff and Lachlan about this.
  
  if (OutputRasters == "Yes"){
    #cl<-makeCluster(Ncores,outfile="")
    #registerDoParallel(cl, cores=Ncores)
    #foreach(i=iter(1:LTR),.packages =c("raster","filematrix") )%dopar%{
    
    for (i in 1:LTR){
      j=TimeNames[i]
      outR<-r
      values(outR)<-TFI_VAL_RASTER_VALS[,i]
      writeRaster(outR,file.path(ResultsDir,"TFI_Rasters",paste0("TFI_VAL",j,".tif")),overwrite=T,datatype="INT1U")
      rm(outR)
      gc()
    }
    #stopImplicitCluster()
    
  }
  return(TFI_Summary)
}











# calculate date sequence BBTFI (accomodating  Hi and Lo fire intesity of first burn to determine TFI) ----------------
#returns list containing  [[1]] the date sequence matrix for each cell of the raster
#                         [[2]] the EFG TFI Lookup for each cell of the raster
#                         [[3]] the raster resolution used.
# throws an error if raster resolution does not match selected FH_ID file

calc_BBTFI<-function(FHanalysis,#the slected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                     cropRasters,
                     TFI_LUT_DF = TFI_LUT#the dataframe read from a csv that gives the lookup table from EVD to MinTFI_LO, MinTFI_HI and MaxTFI
                     ) {
  
  r<-FHanalysis$FH_IDr
  Hectares<-(as.numeric(FHanalysis$RasterRes)/100)^2

  FH_ID<-as.data.frame(values(r))
  names(FH_ID)<-"ID"
  DELWP<-cropRasters$DELWP
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
  
  TYPE_HI<-TYPE==2# TRUE where thetype of the first fire is HI(2)
  TYPE_LO<-TYPE==1# TRUE where thetype of the first fire is LO(1)
  
  BB_LO_TFI_SEASON<-SEAS*(INTS<TFI$MIN_LO_TFI)#if any fire interval is less than the MIN_LO_TFI then the veg has been burnt below TFI multiplication of the true or false by the season returns the season where true and zeor where false 
  BB_LO_TFI_SEASON[BB_LO_TFI_SEASON==0]<-NA#remove zeros because later going to want to calc minimum dates
  BB_HI_TFI_SEASON<-SEAS*TYPE_HI*(INTS<TFI$MIN_HI_TFI)#only in cases where the first firetype is HI and the interval is below MIN_HI_TFI is the veg burnt below TFI
  BB_HI_TFI_SEASON[BB_HI_TFI_SEASON==0]<-NA
  
  #next two line combine dates for fires BBTFI that are below the low TFI threhsold,
  #and those below the Hi TFI theshold( as determined by the TYPE of the first fire)
  BBTFI_COMB<-BB_LO_TFI_SEASON
  BBTFI_COMB[is.na(BB_LO_TFI_SEASON)]<-BB_HI_TFI_SEASON[is.na(BB_LO_TFI_SEASON)]
 
  #separating whether the fire in the second season of the BBTFI was high or low for reporting
  ID<-1:nrow(BBTFI_COMB)
   
    

  BBTFI_ID<-as_tibble(cbind(ID,BBTFI_COMB))%>%gather(Garbage,SEASON,-ID)
  TYPE2_ID<-as_tibble(cbind(ID,TYPE2))%>%gather(Garbage,TYPE,-ID)
 
  ID_LU<-as_tibble(cbind(ID,EFG,FIREFMZ,FIRE_REG,DELWP,PLM))

  BBTFI_BY_TYPE<-right_join(ID_LU,na.omit(cbind(BBTFI_ID[,c("ID","SEASON")],TYPE2_ID[,"TYPE"])))%>%
    count(EFG,FIRE_REG,DELWP,FIREFMZ,PLM,SEASON,TYPE)
  BBTFI_BY_TYPE$Hectares<-BBTFI_BY_TYPE$n*Hectares
  BBTFI_BY_TYPE<-right_join(FIRETYPE_LUT,BBTFI_BY_TYPE)
  BBTFI_BY_TYPE<-right_join(TFI_LUT[,c("EFG","EFG_NAME")],BBTFI_BY_TYPE)
  BBTFI_BY_TYPE<-right_join(FIREFMZ_LUT,BBTFI_BY_TYPE)
  BBTFI_BY_TYPE<-right_join(REG_LUT,BBTFI_BY_TYPE)
  BBTFI_BY_TYPE<-right_join(DELWP_LUT,BBTFI_BY_TYPE)
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

#
makeGS_Sum<-function(TimeSpan = FHanalysis$TimeSpan,
                     writeGSRasters,
                     myLU = GS_LU,
                     myResultsDir = ResultsDir,
                     myCropDetails = cropRasters,
                     myFHResults = FHanalysis,
                     myYSF_LFT = tsf_ysf_mat,
                     writeYears=NULL) {
  PLM<-cropRasters$PLM
  EFG<-cropRasters$EFG
  FIRE_REG<-cropRasters$RGN
  FIREFMZ<-cropRasters$FIREFMZ
  DELWP<-cropRasters$DELWP
  GSYrArray<-NULL
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
    OutTif <-file.path(myResultsDir,"GS_Rasters", paste0("GS_YR_", year, ".tif"))
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
    GSYrArray<-cbind(GSYrArray,Out)
    GSYrres <- c(year, table(Out))
    GSYearSumm <- rbind(GSYearSumm, GSYrres)
    
    print(paste("Finished",year))
    
  }
  colnames(GSYrArray)<-TimeSpan
  myTab<-cbind(EFG,FIRE_REG,DELWP,FIREFMZ,PLM,GSYrArray)
  myTab<-as_tibble(myTab)
  #tablong<-gather(myTab,SEASON,GS,-EFG,-FIRE_REG,-FIREFMZ,-PLM)
  GS_Summary<-myTab%>%
    gather(SEASON,GS,-EFG,-FIRE_REG,-FIREFMZ,-PLM)%>%
    count(EFG,FIRE_REG,DELWP,FIREFMZ,PLM,SEASON,GS)
  rm(myTab)
  gc()
  GS_Summary$Hectares<-GS_Summary$n*cellsToHectares()
  GS_Summary<-left_join(GS_Summary,GS_LUT)
  GS_Summary<-left_join(GS_Summary,TFI_LUT[,c("EFG","EFG_NAME")])
  GS_Summary<-left_join(GS_Summary,FIREFMZ_LUT)
  GS_Summary<-left_join(GS_Summary,REG_LUT)
  return(GS_Summary)
}

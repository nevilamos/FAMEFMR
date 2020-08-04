# calculate date sequence BBTFI (accomodating  Hi and Lo fire intesity of first burn to determine TFI) ----------------
#returns list containing  [[1]] the date sequence matrix for each cell of the raster
#                         [[2]] the EFG TFI Lookup for each cell of the raster
#                         [[3]] the raster resolution used.
# throws an error if raster resolution does not match selected FH_ID file

calcBBTFI_2<-function(FHanalysis,
                      U_AllCombs_TFI=myAllCombs$U_AllCombs_TFI,
                      Index_AllCombs=myAllCombs$Index_AllCombs,
                      TFI_LUT,
                      makeBBTFIrasters=makeBBTFIrasters )#the dataframe read from a csv that gives the lookup table from EVD to MinTFI_LO, MinTFI_HI and MaxTFI)
{
  r<-FHanalysis$FH_IDr
  
  FH_ID<-values(r)
  r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
  
  
  #read the corresponding shapefile for the FireHat analysis chosen
  #convert to a data.frame
  
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  
  INTFields<-names(OutTab)[grep("^INT",names(OutTab))]#fields from FH analysis giving inter fire intervalsof 1:nth fire
  SEASFields<-names(OutTab)[grep("^SEAS",names(OutTab))]#fields giving season of 1:nth fire
  TYPEFields<-names(OutTab)[grep("^FireType",names(OutTab))]#fields giving type of 1:nth fire
  # x<-left_join(FH_ID,OutTab)
  # x<-as.matrix(x)
  # mode(x)<-"integer"
  INTS<-as.matrix(OutTab[,INTFields])#The fire intervals in numbered sequence
  SEAS<-as.matrix(OutTab[,SEASFields[-1]])#The Season of the second fire of each interval to get the season when burning under TFI occurs
  SEAS[SEAS==0]<-NA
  TYPE<-as.matrix(OutTab[,TYPEFields[-1*length(TYPEFields)]])#The type of the first fire of the interval to determine whether high or lo TFI applies
  
  
  TYPE2<-as.matrix(OutTab[,TYPEFields[-1]])#The firetype of the second fire - that determines whether the fire at the date that the second burn event occured was high or low for reporting
  TYPE2[TYPE2==0]<-NA
  
  
  
  TYPE_HI<-TYPE==2# TRUE where thetype of the first fire is HI(2)
  TYPE_LO<-TYPE==1# TRUE where thetype of the first fire is LO(1)
  
  SEAS<-SEAS[U_AllCombs_TFI$FH_ID,]
  INTS<-INTS[U_AllCombs_TFI$FH_ID,]
  TYPE_HI<-TYPE_HI[U_AllCombs_TFI$FH_ID,]
  #TYPE_HI<-TYPE_HI[U_AllCombs_TFI$FH_ID,]
  TYPE1<-TYPE[U_AllCombs_TFI$FH_ID,]
  TYPE2<-TYPE2[U_AllCombs_TFI$FH_ID,]
  
  
  BB_LO_TFI_SEASON<-SEAS*(INTS<U_AllCombs_TFI$MIN_LO_TFI)#if any fire interval is less than the MIN_LO_TFI then the veg has been burnt below TFI multiplication of the true or false by the season returns the season where true and zero where false 
  BB_LO_TFI_SEASON[BB_LO_TFI_SEASON==0]<-NA#remove zeros because later going to want to calc minimum dates
  BB_HI_TFI_SEASON<-SEAS*TYPE_HI*(INTS<U_AllCombs_TFI$MIN_HI_TFI)#only in cases where the first firetype is HI and the interval is below MIN_HI_TFI is the veg burnt below TFI
  BB_HI_TFI_SEASON[BB_HI_TFI_SEASON==0]<-NA
  
  #next two line combine dates for fires BBTFI that are below the low TFI threhsold,
  #and those below the Hi TFI theshold( as determined by the TYPE of the first fire)
  BBTFI_COMB<-BB_LO_TFI_SEASON
  BBTFI_COMB[is.na(BB_LO_TFI_SEASON)]<-BB_HI_TFI_SEASON[is.na(BB_LO_TFI_SEASON)]
  
  #cumulative number of times bbtfi events in each comination of EFG with fire history and grouping polygons
  cumBBTFI<-t(apply(!is.na(BBTFI_COMB),1,FUN = "cumsum"))
  totalTimesBBTFI<-rowSums(!is.na(BBTFI_COMB))
  #totalTimesBBTFI[totalTimesBBTFI==0]<-NA
  #mininimum date in row determines first time that was burned below TFI
  FirstBBTFI<-apply(BBTFI_COMB,1,min,na.rm=T )
  FirstBBTFI[is.infinite(FirstBBTFI)]<-NA
  #need to set cells NA where NA in BBTFI
  cumBBTFI[is.na(BBTFI_COMB)]<-NA
  colnames(cumBBTFI)<-gsub("SEAS","TimesBBTFI",colnames(cumBBTFI))
  
  # have called this combined item BBTFI_WIDEish because is is the same number of rows as the true wide item but cumBBTFI is in sequece for still not spread by year.
  #BBTFI_WIDEish<-cbind(U_AllCombs_TFI,BBTFI_COMB,cumBBTFI,TYPE2,FirstBBTFI,totalTimesBBTFI)
  BBTFI_WIDEish<-cbind(U_AllCombs_TFI,SEAS,cumBBTFI,TYPE2,FirstBBTFI,totalTimesBBTFI)
  BBTFI_LONG<-BBTFI_WIDEish%>%#select(-FirstBBTFI,-totalTimesBBTFI)%>%
    pivot_longer(all_of(names(BBTFI_WIDEish)[(ncol(U_AllCombs_TFI)+1):(ncol(BBTFI_WIDEish)-2)]),
                 names_to = c(".value", "SEQ"),
                 names_pattern = "([^0-9]+)([0-9]+)",
                 values_drop_na=T)
                 
  BBTFI_WIDE<-BBTFI_LONG%>%
    select(-SEQ,-totalTimesBBTFI,-FirstBBTFI)%>%
    pivot_wider(names_from = SEAS,values_from=c(TimesBBTFI,FireType),names_sort=T)
  
  BBTFI_LONG_Summary<- BBTFI_LONG%>%
    group_by(PLM,FIRE_REGION_NAME,DELWP_REGION,EFG_NAME,FIRE_FMZ_NAME,SEAS,FireType,TimesBBTFI)%>%
    summarize(Hectares=sum(Hectares))
  
  if(makeBBTFIrasters){
        values(r)<-totalTimesBBTFI[Index_AllCombs]
        writeRaster(r,file.path(ResultsDir,"TotalTimesBBTFI.tif"),overwrite=T,datatype="INT1U")
        r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
        values(r)<-FirstBBTFI[Index_AllCombs]
        writeRaster(r,file.path(ResultsDir,"FirstBBTFI.tif"),overwrite=T,datatype="INT4U")
        r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
  }
  
  return(list("BBTFI_WIDE"=BBTFI_WIDE,"BBTFI_LONG"=BBTFI_LONG_Summary,"BBTFI_WIDEish"=BBTFI_WIDEish))
  
  
}            
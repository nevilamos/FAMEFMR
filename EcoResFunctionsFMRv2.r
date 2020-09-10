
#Functions used in EcoRes1 calculation of Fire History Interval and spatial realtive abundance related to fire history.
# written by nevil.amos@delwp.vic.gov.au


################################################################################
#revove all empty  directories in a path----------------------------

removeEmptyDirs<-function(rootDir="./Results"){
  
  dirList <- list.dirs(rootDir)[-1]#makes sure the root results direcotry is not deleted
  if(length(dirList)>0){
    delList<-dirList[sapply(dirList, function(x) length(list.files(x))==0)]
    while((length(delList)>0) & (length(dirList)>0)){ 
      unlink(delList,recursive=T)
      dirList<-list.dirs(rootDir)[-1]
      delList<-dirList[sapply(dirList, function(x) length(list.files(x))==0)]
    }
  }
}

# add xy string -----------------------------------------------------------


#adds string of xy of centroids of polygons to allow for rearrangement of  associated data in flattened fire history
add_xystring<-function(myDF){
  #function to add xystring of polygon centroids to SF polygon object
  require(dplyr)
  Coords<-as.data.frame(do.call(rbind,st_geometry(st_centroid(myDF))))
  names(Coords)<-c("X","Y")
  XYString<-paste(Coords$X,Coords$Y,sep="")
  x<-mutate(myDF,XYString)
  return(x)
  
}

###get_Spp_No#############################################################################
#extracts four or five digit species numbers from  HDM paths
get_Spp_No<-function(x="Vector of Sp file Pathnames"){
  Fnames= basename(x)
  pos = regexpr("[0-9][0-9][0-9][0-9]+", Fnames)
  myEnd=pos-1+attr(pos,"match.length")
  y=as.numeric(substr(Fnames,pos,myEnd))
  return(y)
}


###makeCropDetails############################################################################
# function to get the minimum bounding box of the cells with non NA values in a raster and save them to crop other rasters to same extent
#also creates some rasters cropped to correct extent for instance for region and EFG
#also gets indeces of cells to "clip raster of same extent as crop to the shape provided 

makeCropDetails<-function(REG_NO=7,#REG_NO of defined region from input (1:6) or 0 for statewide or 7 for Ad Hoc Poly),
                          RasterRes=RasterRes,
                          PUBLIC_LAND_ONLY="YES",
                          myPoly=clipPoly,#shapefile ofLF_REGIONs( default)or  adhoc region,
                          generalRasterDir = "./InputGeneralRasters"
){
  inputR<-inputRasters(RasterRes)
  inR<-raster(file.path(generalRasterDir,inputR$REGION.tif))
  Template<-inR
  values(Template)<-NA
  #determines which file to use for masking to regions
  if(REG_NO%in%1:6){
    Shape<-read_sf(myPoly)
    Shape<-Shape[Shape$REGION_NO==REG_NO,]
    cn<-cellnumbers(Template,Shape)
    cn<-cn$cell_
    RGN<-Template
    values(RGN)[cn]<-REG_NO
    
  }
  if(REG_NO==99){
    Shape<-read_sf(myPoly)
    cn<-cellnumbers(Template,Shape)
    cn<-cn$cell_
    RGN<-inR
    #values(RGN)[cn]<-REG_NO
  }
  if(REG_NO ==7){
    Shape<-read_sf(myPoly)
    cn<-cellnumbers(Template,Shape)
    cn<-cn$cell_
    RGN<-Template
    values(RGN)[cn]<-REG_NO
  }
  
  
  
  x=RGN
  x.matrix <- is.na(as.matrix(x))
  colNotNA <- which(colSums(x.matrix) != nrow(x))
  rowNotNA <- which(rowSums(x.matrix) != ncol(x))
  Extent <- extent(x,
                   r1 = rowNotNA[1],
                   r2 = rowNotNA[length(rowNotNA)],
                   c1 = colNotNA[1],
                   c2 = colNotNA[length(colNotNA)])
  
  
  RGN_ras<-crop(RGN,Extent)
  FIREFMZ_ras<-crop(raster(file.path(generalRasterDir,inputR$FIREFMZ.tif)),Extent)
  DELWP_ras<-crop(raster(file.path(generalRasterDir,inputR$DELWP.tif)),Extent)
  EFG_ras<-crop(raster(file.path(generalRasterDir,inputR$EFG.tif)),Extent)
  IDX<-values(crop(raster(file.path(generalRasterDir,inputR$IDX.tif)),Extent))
  PLM_ras<-crop(raster(file.path(generalRasterDir,inputR$PLM_GEN.tif)),Extent)
  EFG_ras<-mask(EFG_ras,RGN_ras)
  # if choice has been made to restrict to public land ( default) then EFG is masked to public land
  if(PUBLIC_LAND_ONLY == "YES"){
    EFG_ras<-mask(EFG_ras,PLM_ras)
    DELWP_ras<-mask(DELWP_ras,PLM_ras)
    FIREFMZ_ras<-mask(FIREFMZ_ras,PLM_ras)
  }
  PLMVals<-as.integer(values(PLM_ras))
  EFGvals<-as.integer(values(EFG_ras))
  RGNvals<-as.integer(values(RGN_ras))
  FIREFMZvals<-as.integer(values(FIREFMZ_ras))
  DELWPvals<-as.integer(values(DELWP_ras))
  
  
  CropDetails<-list("Raster" = RGN_ras,
                    "Extent"=Extent,
                    "clipIDX"=cn,
                    "EFG"=EFGvals,
                    "RGN"=RGNvals,
                    "IDX"=IDX,
                    "DELWP"=DELWPvals,
                    "FIREFMZ"=FIREFMZvals,
                    "PLM"=PLMVals)
  return(CropDetails)
  
}





###notAllIN#######
#function to check whether all values in vector are in in another vector of permitted values, retuns fals
notAllIn<-function(x,v=V){
  anyNA(unlist(lapply(x,match,v)))
}
###FHProcess###############################################################################################################################

#function converts rawFH dataset into flattend format ( single non overlapping polygon attributed with unique fire sequences) calculates all inter-fire intervals time since fire for a sequence of years and output as rasters
FHProcess<-function(rawFH ="path of the rawFH file to use - a shapefile",
                    start.SEASON=NULL, #"first season for which output is wanted ( four digit year as integer)" if NUll then second season in in history  is used cannot use first season because it has no interval, this may still fail if there is no overlap,
                    end.SEASON=NULL, #"last season required if NULL then largest value in fire history scenario used"
                    OtherAndUnknown =  2# (2,1,NA)value to use for cases where fire type is "OTHER" or "UNKNOWN",1 ="BURN",2="BUSHFIRE",NA = Fire excluded from analysis default is 2 ("BUSHFIRE")
                    
){
  myDF<-st_read(rawFH)
  # check that the input shapefile (rawFH) contains the two required fields and that these do not have a missing values
  myDFNames<-names(myDF)
  if((!"SEASON"%in% myDFNames)) stop ("raw FH does not contain field 'SEASON'")
  if(anyNA(myDF$SEASON)) stop ("rawFH has missing values in  field 'SEASON'")
  if((!"FIRETYPE"%in% myDFNames)) stop ("raw FH does not contain field 'SEASON'")
  if(notAllIn(x=myDF$FIRETYPE,v=validFIRETYPE)) stop ("rawFH has missing  or invalid values in  field 'FIRETYPE'")
  
  #timespan ( range of consecutive years) for which Fire History sequences are calculated
  #start and end season  for calcuation is determined by the min+1 and max values in the dataset,
  #or if not null from start and end seasons defined in the settings file.
  # the earliest start season is the date of a second season in the input data, next 3 lines prevent manual start
  #season setting  less than this value, which would cause an error.
  min.SEASON<-sort(unique(myDF$SEASON))[2]# second season in fire history
  if(is.null(start.SEASON)){
    start.SEASON=min.SEASON
  } else {
    if(start.SEASON<min.SEASON){
      start.SEASON=min.SEASON
    }else{
      start.SEASON=start.SEASON
    }
    
  }
  
  if(is.null(end.SEASON)){
    max.SEASON<-max(myDF$SEASON)
  } else {
    max.SEASON = end.SEASON
  }
  
  TimeSpan<-start.SEASON:max.SEASON
  
  #recodes FIRETYPE  to integer values value for Other and Unkonwn is set in settings file
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="BURN"]<-1
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="BUSHFIRE"]<-2
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="OTHER"]<-OtherAndUnknown
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="UNKNOWN"]<-OtherAndUnknown
  
  #adds a field containing a string representation of the x and y values of the
  #centroid of each polygon in myDF to use as a unique identifier for unique
  #polygons in the FH analysis ( ie two sptailly identical polygons with
  #different SEASON of Firetype values will share an uniqe xycoords string) this
  #is used to get the wide format sequence for each spatially unique polygon )
  myDF<-add_xystring(myDF)
  
  myDF<-unique(myDF[,c("geometry",
                       "XYString",
                       "SEASON",
                       "FIRETYPE_NO")])
  #orders the spatially unique polygons by SEASON then firetype and adds a
  #sequential number to them , this is what allows the flattening  and subsequent
  #reduction to fire sequences witout interveneing no fire year
 
  myDF<-myDF[with(myDF,order(XYString, SEASON, FIRETYPE_NO)),] 
  
  myDF$Sequence<-1 
  myDF$Sequence<-unlist(lapply(split(myDF$Sequence, myDF$XYString), cumsum)) 
  
  #reduction to the unique sequences of fires with no gaps using dplyr::pivot_wider 
  # in order to acheive this efficently the sf object made from the rawFH file is converted to a standard data frame iwtha geometry column then coverted back to an sf object at the end of the process
  tic("Making OutDF")
  OutDF<-myDF%>%
    select(XYString,Sequence,FireType=FIRETYPE_NO,SEAS=SEASON)%>%
    mutate(Sequence=sprintf("%02d",Sequence))%>%
    as.data.frame%>%
    pivot_wider(names_from = Sequence,values_from=c(FireType,SEAS),names_sep="")
  # get just the names for the season  and FireType sequence columns
  SEASNames<-names(OutDF)[grep(pattern = "SEAS",names(OutDF))]
  FTNames<-names(OutDF)[grep(pattern = "FireType",names(OutDF))]
  #Calculates all inter-fire intervals note by offsetting matrices by one column
  SEAS_Matrix<-as.matrix(OutDF[,SEASNames])
  FT_matrix<-as.matrix(OutDF[,FTNames])
  Cols<-ncol(SEAS_Matrix)
  SEAS_Matrix[SEAS_Matrix==0]<-NA
  Interval<-SEAS_Matrix[,2:Cols]-SEAS_Matrix[,1:Cols-1]
  IntNames<-paste("INT",sprintf("%02d",1:(Cols-1)),sep="")
  colnames(Interval)<-IntNames
  
  #Binds the intevals calcuated back to the main dataframe
  OutDF<-cbind(OutDF,Interval)%>%
    select(-XYString)
  
  OutDF<-st_as_sf(OutDF)
  toc()
 
  #calcuating the last burnt season for each sequences for each year 
  #this process is duplicated here and in CALC_TFI2 function - ideally it should only be run once but it is pretty fast so probably does not matter tooo much.
  #
  LTR<-length(TimeSpan)
  SEAS_Matrix[is.na(SEAS_Matrix)]<-0
  LBY<-matrix(NA,nrow(SEAS_Matrix),LTR)
  for(i in 1:LTR){
    try({
      y=TimeSpan[i]
      LBY[,i]<-LBY_f(M=SEAS_Matrix,y)
      #print(y)
      
    })
  }
  
  #subtraction leads to matrix transpostion which needs to be retransposed in the next couple of lines
  tYSF<-TimeSpan-t(LBY)
  YSF<-t(tYSF)
  
  YSFNames<-paste0("YSF",TimeSpan)
  LBYNames<-paste0("LBY",TimeSpan)
  LFTNames<-paste0("LFT",TimeSpan)
  colnames(YSF)<-YSFNames
  colnames(LBY)<-LBYNames
  
  
  tic("calculating lookup matrix for getting last firetype by year")
  SEAS_Matrix[SEAS_Matrix==0]<-NA
  LUM<-matrix(NA,nrow(SEAS_Matrix),max.SEASON)
  for (i in 1:nrow(SEAS_Matrix)){
    R<-i
    C<-as.numeric(na.omit(SEAS_Matrix[i,]))
    V<-(FT_matrix[i,(1:length(C))])
    LUM[R,C]<-V
  }
  toc()
  
  tic("calculating last fire type")
  LFT<-matrix(NA,nrow(SEAS_Matrix),LTR)
  for(i in 1:nrow(SEAS_Matrix)){
    LFT[i,]<-LUM[i,LBY[i,]]}
  colnames(LFT)<-LFTNames
  toc()
  
  OutDF<-cbind(OutDF,YSF)
  OutDF<-cbind(OutDF,LBY)
  OutDF<-cbind(OutDF,LFT)
  OutDF$ID<-as.integer(rownames(OutDF))
  print("completed making FH object")
  
  
  
  
  
  gc()
  results<-list(OutDF=OutDF,TimeSpan=TimeSpan,  YSFNames=YSFNames,LBYNames=LBYNames, LFTNames=LFTNames)
  return(results)
}




###makeLU_List#############################################################################################
# function creates a list of Lookup arrays for each taxon(VBA_CODE/TAXON_ID) for  YSF x EFGNO x FireTypeNo these are then used in
makeLU_List<- function(myHDMSpp_NO=HDMSpp_NO,#
                       myAbundDataLong=ExpertDataLong){
  
  
  
  myList <- list() 
  for(i in myHDMSpp_NO){
    
    y<-myAbundDataLong[myAbundDataLong$VBA_CODE==i,]
    b=(y$YSF)+1
    c=y$EFG_NO
    d=y$FireTypeNo
    e=y$Abund
    x<-array(NA,dim=c(max(b),40,4))
    for(j in 1:nrow(y)){
      x[b[j],c[j],d[j]]<-e[j]
    }
    myList[[as.character(i)]]<-x
    rm(y)
  }
  return(myList)
}

###makeSppYearSum2#######################################################
#function returns summary of species summed relative abundances by year
#also if writeSpRasters==TRUE it writes species raster for each year (slows processing)
#may require excessive ram to run statewide at 75m resolution.
makeSppYearSum2<-function(FHanalysis,
                          myHDMSpp_NO = HDMSpp_NO,
                          writeSpRasters = writeSpRasters,
                          myLU_List = LU_List,
                          ResultsDir = ResultsDir,
                          HDMVals = HDMVals,
                          TaxonList = TaxonList,
                          writeYears=NULL,
                          writeSp =NULL
) {
  
  TimeRange<-as.integer(FHanalysis$TimeSpan)
  TimeNames<-as.character(FHanalysis$TimeSpan)
  LTR<-length(TimeRange)
  r<-FHanalysis$FH_IDr
  
  r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
  TimeSpan<-FHanalysis$TimeSpan
  myDF<-FHanalysis$OutDF
  #remove geometry FHanalysis DF to create a standard dataframe so columns can be subset without "stciky geometry of origina spatial Features data frame
  st_geometry(myDF)<-NULL
  
  #create empty matrix with rownames containing HDM TAXON ID numbers and column names of the years ( SEASONS) to house the Species abundance data for each year.
  SpYearSumm<-matrix(NA,nrow=(length(myHDMSpp_NO)),ncol=LTR,dimnames=list(as.character(myHDMSpp_NO),TimeNames))
  #loop through calcualtion of per cell speices abundance values and if flagged output of rasters of these values.
  for (sp in myHDMSpp_NO) {
    print(sp)
    mySpp<-as.character(sp)
    #get the lookup array of abundance values from the list of species value lookup 
    LU = myLU_List[[mySpp]]
    #gets the HDMVals for the relevant species from the matrix of HDMvalues by species. values multiplied by 100 so that they can later be converted to integer if necessarywithout losing small values
    HDM_Vector<-as.vector(HDMVals[,grep(as.character(sp),colnames(HDMVals))])*100
    #makes matrices of YSF,EFG,and LFT to use in lookup from LU array
    YSF_M<-as.matrix(myDF[myAllCombs$U_AllCombs_TFI$FH_ID,FHanalysis$YSFNames])+1
    LFT_M<-as.matrix(myDF[myAllCombs$U_AllCombs_TFI$FH_ID,FHanalysis$LFTNames])
    EFG_M<-matrix(myAllCombs$U_AllCombs_TFI$EFG,nrow(YSF_M),ncol(YSF_M))
    #looks up the cell-wise psecies values for for the species abundance balues by ineces in array
    Spp_M<-array(LU[cbind(as.vector(YSF_M),as.vector(EFG_M),as.vector( LFT_M))],dim=dim(YSF_M))
    #Multiplies these values by cell-wise HDM values - effectively masking out values where the species does not occur.
    Spp_Val_Cell_Year<-Spp_M[myAllCombs$Index_AllCombs,]*HDM_Vector
    colnames(Spp_Val_Cell_Year)<-TimeNames
    #get the sum of cell values for each year for the species
    #put them in the compilation data frame
    SpYearSumm[mySpp,]<-colSums(Spp_Val_Cell_Year,na.rm=T)
    
    gc()
    #if writing species rasters is desired ( set in settings)
    #next section writesout species rasters
    #this is by far the most time consuming part of the FAME processing
    
    if(writeSpRasters=="Yes"){
      for(myYear in as.character(writeYears))
        if (sp%in%writeSp|is.null(writeSp)){
          emptySpraster <- r
          values(emptySpraster) <- Spp_Val_Cell_Year[,myYear]
          writeRaster(emptySpraster,
                      OutTif,
                      options = c("COMPRESS=LZW", "TFW=YES"),
                      datatype = 'INT1U',
                      overwrite = TRUE
          )
        }
    }
    
    
  }
  #join species details from TaxonList  to the output tables
  SpYearSumm<-rownames_to_column(as.data.frame(SpYearSumm))
  names(SpYearSumm)[1]<-"TAXON_ID"
  TL<-TaxonList%>%
    mutate(TAXON_ID=as.character(TAXON_ID))
  SpYearSummWide<-right_join(TL,SpYearSumm)
  write.csv(SpYearSummWide,file.path(ResultsDir,"SpYearSummWide.csv"))
  SpYearSummLong<-right_join(TL,SpYearSumm%>%pivot_longer(-TAXON_ID,names_to="SEASON",values_to="SUM_RAx100"))
  write.csv(SpYearSummLong,file.path(ResultsDir,"SpYearSummLong.csv"))
  
  return(SpYearSummWide)
  }


###calcDeltaAbund#######################################################################################

#Calculates baseline RA based on input baseline years and deviation from  a
#baseline ( either a single year or a several years, usually a sequence eg
#1980-1990) for each  year, output written to CSV file.
#"SppSummChangeRelativetoBaseline.csv" 

calcDeltaAbund<- function(SpYearSummSpreadbyYear =SpYearSummWide,
                          TimeSpan =FHanalysis$TimeSpan,
                          myBaseline,
                          ResultsDir,
                          HDMSpp_NO,
                          TaxonList)
{

  #to get % of baseline need to define which columns provide the baseline ( one or mean of several using apply (mean)) then divide remaining values by this column.
  if (length(myBaseline==1)){
    Baseline<-SpYearSummSpreadbyYear[,as.character(myBaseline)]
  } else {
    Baseline<-apply(SpYearSummSpreadbyYear[,as.character(myBaseline)],1,mean)
  }
  SpYearSummSpreadbyYear$Baseline<-Baseline
  # gets the integer value for current year used in next line so that changes to baseline
  #are only displayed for future years or if no future years then years since baseline.
  ThisYear<-as.integer(format(Sys.Date(), "%Y"))
  SinceYear<-ifelse(sum(TimeSpan>ThisYear)>0,
                    ThisYear,
                    max(myBaseline))
  
  #calcualtes th changes from baseline
  Deltas<-as.matrix(SpYearSummSpreadbyYear[,as.character(TimeSpan[TimeSpan>SinceYear])]/Baseline)
  names(Deltas)<-paste(names(Deltas),"prop baseline")
  
  # calcuates two potential comparivie metrics , the total number of years below threshold and whether the last year is below threshold
  NoLessthanThreshhold<-rowSums(Deltas<=SpYearSummSpreadbyYear$CombThreshold)
  LastLessThanThreshold<-Deltas[,ncol(Deltas)]<=SpYearSummSpreadbyYear$CombThreshold
  
  #Subsets input columns and appends results to make an output table
  ChangeRelativeToBaseline<-cbind(SpYearSummSpreadbyYear[,c("TAXON_ID",
                                                            #"ShortName",
                                                            "COMMON_NAME",
                                                            "NAME",
                                                            "DIVNAME",
                                                            "EPBC_ACT_STATUS",
                                                            "VIC_ADVISORY_STATUS",
                                                            "CombThreshold")],
                                  Deltas,
                                  NoLessthanThreshhold,
                                  LastLessThanThreshold)
  #wirtes the table to file
  write.csv(ChangeRelativeToBaseline,
            file.path(ResultsDir,"SppSummChangeRelativetoBaseline.csv"),
            row.names=F)
  return(ChangeRelativeToBaseline)
  
}




###calcDraftSpList#####################################################################################################

#Calculate the proportion of cells for the HDM in the region for each species works by using the indices of the standard dimesions 
#raster that are in the supplied shapefile region boundary
calcDraftSpList<-function(REG_NO,#REG_NO of defined region from input (1:6) or 0 for statewide or 7 for Ad Hoc Poly),
                          RasterRes=225,
                          PUBLIC_LAND_ONLY="YES",
                          myPoly=myPoly,#shapefile ofLF_REGIONs( default)or  adhoc region,
                          generalRasterDir = "./InputGeneralRasters",
                          splist ="./ReferenceTables/DraftTaxonListStatewidev2.csv",
                          HDMVals=HDMVals225){
  load(HDMVals)
  REG_NO<-as.integer(as.numeric(REG_NO))
  splist<-read.csv(splist)
  CropDetails<-makeCropDetails(REG_NO=REG_NO,
                               RasterRes=RasterRes,
                               PUBLIC_LAND_ONLY=PUBLIC_LAND_ONLY,
                               myPoly=myPoly,
                               generalRasterDir = "./InputGeneralRasters"
  )
  
  cellsInArea<-colSums(HDMVals[CropDetails$clipIDX,])
  cellsInState<-colSums(HDMVals)
  areaProp<-signif(cellsInArea/cellsInState,digits = 2)
  TAXON_ID<-as.numeric(colnames(HDMVals))
  myDF<-data.frame(TAXON_ID,cellsInState,cellsInArea,areaProp)
  myDF<-left_join(splist,myDF)
  return(myDF)
}

###calcSppEFGLM#####################################################################################

#Calculate the species in each EFG in given area for GSO calcuations. works by using the indices of the standard dimesions 
#raster that are in the supplied shapefile region boundary ( via function makeCropDetails)
calcSpp_EFG_LMU<-function(REG_NO,#REG_NO of defined region from input (1:6) or 0 for statewide or 7 for Ad Hoc Poly),
                          RasterRes=225,
                          PUBLIC_LAND_ONLY="YES",
                          myPoly=myPoly,#shapefile ofLF_REGIONs( default)or  adhoc region,
                          generalRasterDir = "./InputGeneralRasters",
                          splist ="./ReferenceTables/DraftTaxonListStatewidev2.csv",
                          HDMVals=HDMVals225,
                          EFGRas=EFGRas,
                          TFI_LUT=TFI_LUT){
  options(stringsAsFactors = F)
  load(HDMVals)
  mySpList<-read.csv(splist)[,c( "TAXON_ID","COMMON_NAME","NAME")]
  EFG<-values(raster(EFGRas))
  REG_NO<-as.integer(as.numeric(REG_NO))
  CropDetails<-makeCropDetails(REG_NO=REG_NO,
                               RasterRes=RasterRes,
                               PUBLIC_LAND_ONLY=PUBLIC_LAND_ONLY,
                               myPoly=myPoly,
                               generalRasterDir = "./InputGeneralRasters"
  )
  
  EFG<-EFG[CropDetails$clipIDX]
  EFG[is.na(EFG)]<-99
  HDMVals<-HDMVals[CropDetails$clipIDX,]
  mode(EFG)<-"integer"
  A<-aggregate.Matrix(HDMVals,EFG,fun='sum')
  myDf<-as.data.frame(as.matrix(A))
  myDf$EFG_NO<-as.integer(rownames(myDf))
  myDf<-gather(myDf,key = "TAXON_ID","CellCount",-EFG_NO)
  myDf<-myDf[myDf$CellCount>0,]
  myDf$TAXON_ID<-as.integer(myDf$TAXON_ID)
  myDf$EFG<-as.integer(myDf$EFG)
  myDf$ha<-myDf$CellCount*5.0625
  myDf<-left_join(myDf,TFI_LUT[,c("EFG","EFG_NAME")],by=c("EFG_NO"="EFG"))
  myDf<-left_join(myDf,mySpList)
  EFG_AREAS<-as.data.frame(table(EFG))
  EFG_AREAS$ha<-EFG_AREAS$Freq*((RasterRes/100)^2)
  EFG_AREAS$EFG<-as.numeric(levels(EFG_AREAS$EFG))
  EFG_AREAS<-right_join(TFI_LUT[,c("EFG","EFG_NAME")],EFG_AREAS,by="EFG")
  write.csv(EFG_AREAS,file.path(ResultsDir,"EFG_AREAS.csv"),row.names = F)
  write.csv(myDf,file.path(ResultsDir,"Spp_EFG_LMU.csv"),row.names = F)
}

###inputRasters-----------------------------------------------------------------------

#get names of input rasters depending on cellSize
inputRasters<-function(x=RasterRes){
  #General Input Rasters change name depending on Raster Res
  if (x==225){
    REGION.tif<-"LF_REGION_225.tif"
    EFG.tif<-"EFG_NUM_225.tif"
    PLM_GEN.tif<-"PLM_GEN_225.tif"
    IDX.tif<-"IndexVals225.tif"
    FIREFMZ.tif<-"FIRE_FMZ_225.tif"
    DELWP.tif<-"DELWP_REGION_225.tif"
    
  }else{  
    REGION.tif<-"LF_REGION_75.tif"
    EFG.tif<-"EFG_NUM_75.tif"
    PLM_GEN.tif<-"PLM_GEN_75.tif"
    IDX.tif<-"IndexVals75.tif"
    FIREFMZ.tif<-"FIRE_FMZ_75.tif"
    DELWP.tif<-"DELWP_REGION_75.tif"
  }
  y <-list("REGION.tif"=REGION.tif,
           "EFG.tif"=EFG.tif,
           "PLM_GEN.tif"=PLM_GEN.tif,
           "IDX.tif"=IDX.tif,
           "FIREFMZ.tif"=FIREFMZ.tif,
           "DELWP.tif"=DELWP.tif)
  return(y)
}

# Function joins Lookup tables to DF containing ID_NO :Name combinations at the moment the LUTS are hard wired.

Join_Names<-function(myDF,#dataframe or similar containing indices for the LUTS listed 
                     LUTS=c("TFI_LUT","TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")){
  for(i in LUTS){
    try(myDF<-left_join(myDF,get(i)))}
    return (myDF)
}

# Function calclulate multiplier form cells to hecatres for cell resolution in Metres ( usually from RasterRes in settings file)----------
cellsToHectares<-function(RasterMetres=RasterRes){
  (RasterMetres/100)^2
}


###FUNCTIONS NOT CURRENTLY USED ####################################################
#st_parallel--------------------------------------------------------
#possible parallel version of st_ functions  not yet expolored from
# https://www.spatialanalytics.co.nz/post/2017/09/11/a-parallel-function-for-spatial-analysis-in-r/ 
#and https://www.spatialanalytics.co.nz/post/2018/04/01/fixing-st-par/
#for linux only
#latest info on this on github suggests that there are not likely to be eay spoeed gains of parallell at the moment (July 2020)
# https://github.com/r-spatial/sf/issues/611

# Paralise any simple features analysis.
st_parallel <- function(sf_df, sf_func, n_cores, ...){
  
  # Create a vector to split the data set up by.
  split_vector <- rep(1:n_cores, each = nrow(sf_df) / n_cores, length.out = nrow(sf_df))
  
  # Perform GIS analysis
  split_results <- split(sf_df, split_vector) %>%
    parallel::mclapply(function(x) sf_func(x, ...), mc.cores = n_cores)
  
  
  # Define the output_class. If length is greater than two, then grab the second variable.
  output_class <- class(split_results[[1]])
  if (length(output_class) == 2){
    output_class <- output_class[2]
  }
  
  # Combine results back together. Method of combining depends on the output from the function.
  if (output_class == "matrix"){
    result <- do.call("rbind", split_results)
    names(result) <- NULL
  } else if (output_class == "sfc") {
    result <- do.call("c", split_results)
    result <- sf_func(result) # do.call combines the list but there are still n_cores of the geometry which had been split up. Running st_union or st_collect gathers them up into one, as is the expected output of these two functions. 
  } else if (output_class %in% c('list', 'sgbp') ){
    result <- do.call("c", split_results)
    names(result) <- NULL
  } else if (output_class == "data.frame" ){
    result <- do.call("rbind", split_results)
  } else {
    stop("Unknown class. st_parallel only accepts the following outputs at present: sfc, list, sf, matrix, sgbp.")
  }
  
  # Return result
  return(result)
}



###makeYSF_LFT_YEAR_RASTERS#########################################################################################
#function to export rasters of TSF and YSF for all years # optional rest of script does not depend on these being made.
#this is not used in FAME v 1 or 1.1

makeYSF_LFT_YEAR_RASTERS<-function(myFHResults = FHResults,
                                   myCropDetails=CropDetails,
                                   myYSF_TSF_Dir =YSF_TSF_Dir){
  print ("Parallel lookup to FH_ID to make rasters of FT and YSF")
  r<-raster(FH_ID.tif)
  
  #reduced the number of cores used here because was running out of ram
  #Ncores<-8
  cl<-makeCluster(Ncores)
  registerDoParallel(cl, cores=Ncores)
  r<-crop(r,myCropDetails$Extent)
  v <- values(r)
  #r1<-myCropDetails$Raster
  #values(r1)<-NA
  OutTab<-myFHResults$OutDF
  st_geometry(OutTab)<-NULL
  
  LU_Names<-c(myFHResults$YSFNames,myFHResults$LFTNames)
  LU<-as.matrix(OutTab[v,LU_Names])
  mode(LU)<-"integer"
  
  
  foreach(Col=iter(LU_Names),.packages = "raster")%dopar%{
    out<-myCropDetails$Raster
    values(out)<-LU[,Col][v]
    writeRaster(out,file.path(myYSF_TSF_Dir,paste0(Col,".tif")),options=c("COMPRESS=LZW", "TFW=YES"),datatype='INT2S', overwrite=TRUE)
    rm(out)
  }
  
  stopCluster(cl)
}


###makeYSF_LFT_matrix#######################################################################
#make matrix of cell wise values of YSF and LFT was used in makeSpYearSumm fucntion in Fame v1 and 1.1
#this is not required by the new makeSpYearSumm2 fucntion that replaces it in FAME V2
makeYSF_LFT_matrix<-function(FHanalysis = FHanalysis,
                             myCropDetails=CropRasters,
                             FH_ID.tif=FHanalysis$FH_IDr){
  r<-FHanalysis$FH_IDr
  r<-crop(r,myCropDetails$Extent)
  v <- values(r)
  #r1<-myCropDetails$Raster
  #values(r1)<-NA
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  
  LU_Names<-c(FHanalysis$YSFNames,FHanalysis$LFTNames)
  LU<-data.matrix(OutTab[v,LU_Names],rownames.force = F)
  mode(LU)<-"integer"
  return(LU)
}

###wLog###############################################################################

# function to write Variable name and variable values to logFile ( specified in script ) to keep record of options chosen during the process.
#wLog<-function(x,y=myLogFile){
#  write(paste(deparse(substitute(x)), "=" ,x),y,append =T)
#}

###makeSppYearSum#######################################################
#function returns summary of species summed relative abundances by year
#also if writeSpRasters==TRUE it writes species raster sor each year (slows processing)
#works as fast as old foreach version without blowing out RAM for statewide 225m
makeSppYearSum<-function(TimeSpan = rv$FHanalysis$TimeSpan,
                         myHDMSpp_NO = HDMSpp_NO,
                         writeSpRasters = writeSpRasters,
                         myLU_List = LU_List,
                         YSF_TSF_Dir = YSF_TSF_Dir,
                         ResultsDir = ResultsDir,
                         EFG = cropRasters$EFG,
                         myCropDetails = cropRasters,
                         HDMVals = HDMVals,
                         myFHResults = FHanalysis,
                         myYSF_LFT = tsf_ysf_mat,
                         TaxonList = myTaxonList,
                         writeYears=NULL,
                         writeSp =NULL) {
  SpYearSumm <- NULL
  for (year in TimeSpan) {
    print(paste("Starting sp abund for",year))
    myYSF <- paste0("YSF", year)
    YSF <- myYSF_LFT[, myYSF] + 1
    myLFT <- paste0("LFT", year)
    LFT <- myYSF_LFT[, myLFT]
    LFT[LFT == 0] <- NA
    myDim <- length(YSF)
    Mask_idx <- (1:myDim)
    RegMaskVal <- YSF + EFG + LFT    #+RGN    
    M<-cbind(YSF,EFG,LFT)
    for (sp in myHDMSpp_NO) {
      #print(sp)
      LU = myLU_List[[as.character(sp)]]
      SpMask <- HDMVals[, as.character(sp)]
      SpMask[SpMask==0]<-NA #this row is needed so that next evaluates the masking cells correctly because NA + FALSE=FALSE not NA
      getVals <- Mask_idx[!is.na(RegMaskVal + SpMask)]#
      OutTif <-file.path(ResultsDir,"RA_Rasters", paste0("Sp_", sp, "_YR_", year, ".tif"))
      #OutName<-paste0(year,"_",sp)
      
      Out <- array(NA, myDim)
      
      
      Out[getVals]<-as.integer(LU[M[getVals,]]*100)
      
      if(writeSpRasters=="Yes"){
        
        
        if(year%in%writeYears|is.null(writeYears)){
          if (sp%in%writeSp|is.null(writeSp)){
            emptySpraster <- myCropDetails$Raster
            values(emptySpraster) <- as.vector(Out)
            writeRaster(
              emptySpraster,
              OutTif,
              options = c("COMPRESS=LZW", "TFW=YES"),
              datatype = 'INT1U',
              overwrite = TRUE
            )
          }
        }
      }
      spYrres <- c(sp, year, sum(Out, na.rm = T))
      SpYearSumm <- rbind(SpYearSumm, spYrres)
    }
    print(paste("Finished",year))
    
  }
  colnames(SpYearSumm)<-c("TAXON_ID",	"SEASON",	"SUM_RAx100")
  SpYearSumm<-as.data.frame(SpYearSumm)
  SpTotals<-SpYearSumm%>%group_by(TAXON_ID)%>%summarize(total = sum(SUM_RAx100))
  SpGterThan0<-SpTotals$TAXON_ID[SpTotals$total>0]
  #print(names(SpYearSumm))
  #print(head(SpGterThan0))
  SpYearSumm<-SpYearSumm[SpYearSumm$TAXON_ID%in%SpGterThan0,]
  
  SpYearSumm<-left_join(SpYearSumm,TaxonList)
  return(SpYearSumm)
}

###makeHDMVals######################################################################
#makes a matrix of HDM values(1,NA) constrained to those cells that are indexed in the cropped area 
makeHDMValsfromRasters<-function(myHDMSpp_NO=HDMSpp_NO,
                                 myCropDetails=cropRasters
){
  HDMPaths<-dir(myCropDetails$HDM_RASTER_PATH,full.names=T,pattern =".tif$")
  HDMPaths<-HDMPaths[get_Spp_No(HDMPaths)%in%myHDMSpp_NO]
  
  print("reading HDMvalues")
  
  #Ncores<-8
  
  cl<-makeCluster(Ncores)
  registerDoParallel(cl, cores=Ncores)
  myHDMVals<-foreach(i=iter(HDMPaths),.combine = cbind,.packages="raster")%dopar%{
    myVals<-values(raster(i))[myCropDetails$IDX]
    myVals}
  stopCluster(cl)
  
  colnames(myHDMVals)<-as.character(get_Spp_No(HDMPaths))
  return(myHDMVals)
}
makeHDMVals<-function(myHDMSpp_NO=HDMSpp_NO,
                      myCropDetails=cropRasters,
                      RasterRes=FHanalysis$RasterRes){
  load(paste0("./HDMS/HDMVals",RasterRes,".rdata"))
  myHDMVals<-HDMVals[myCropDetails$IDX,as.character(myHDMSpp_NO)]
  return(myHDMVals)
}

# rawTiffRead -----------------------------------------------------------

#faster reading of values of Tiff file to vector than provided by raster package (needs the nodata value to be entered).
rawTiffRead<-function(x=singlebandTiff,y=tiff_na_value){
  options(warn=-1)
  z<-as.vector(readTIFF(x,as.is = T))
  z[z==y]<-NA
  return(z)
  options(warn=0)
}

###makeSummaryGraphs####################################################################################################
#provides  graphical summary of relative abundance by year( change in specie relative abundance over time)  and a plotly html file that allows hover over to show species names for each line
#Maybe this should be made availble in the shiny ap as a default screen output.  Ideally with ability to select indivdual species lines for display  cannot work out how to define dropdown (species name(s) selected) dynamically)
makeSummaryGraphs<-function(SpYearSumm,
                            ResultsDir,
                            HDMSpp_NO,
                            outputFH){
  
  No_Of_Species<-length(unique(SpYearSumm$TAXON_ID))
  pal=rainbow(10)
  myPlot <- ggplot(data = SpYearSumm,
                   aes(x = SEASON,
                       y = SUM_RAx100,
                       color=COMMON_NAME)) + 
    geom_line()+
    theme(legend.position="none") +
    ggtitle(paste(outputFH,"\n",No_Of_Species,"Species\n"))
  
  
  ggp<-ggplotly(myPlot,tooltip = "COMMON_NAME") 
  htmlwidgets::saveWidget(as_widget(ggp), file.path(getwd(),ResultsDir,"SpYearSummGraph.html"))
  #unlink(file.path(ResultsDir,SpYearSummGraph_files),recursive = T)
}
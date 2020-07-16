
#Functions used in EcoRes1 calculation of Fire History Interval and spatial realtive abundance related to fire history.
# written by nevil.amos@delwp.vic.gov.au


################################################################################
#revove all empty  directories in a path----------------------------

removeEmptyDirs<-function(rootDir="./results"){
  
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
# rawTiffRead -----------------------------------------------------------

#faster reading of values of Tiff file to vector than provided by raster package (needs the nodata value to be entered).
rawTiffRead<-function(x=singlebandTiff,y=tiff_na_value){
  options(warn=-1)
  z<-as.vector(readTIFF(x,as.is = T))
  z[z==y]<-NA
  return(z)
  options(warn=0)
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




###wLog###############################################################################

# function to write Variable name and variable values to logFile ( specified in script ) to keep record of options chosen during the process.
wLog<-function(x,y=myLogFile){
  write(paste(deparse(substitute(x)), "=" ,x),y,append =T)
}
###FHProcess###############################################################################################################################

#function calculates all inter-fire intervals time since fire for a sequence of years and output as rasters
FHProcess<-function(flattenedFH ="path of the FlattenedFH file to use - a shapefile",
                    start.SEASON=NULL, #"first season for which output is wanted ( four digit year as integer)" if NUll then second season in in history  is used cannot use first season because it has no interval, this may still fail if there is no overlap,
                    end.SEASON=NULL, #"last season required if NULL then largest value in fire history scenario used"
                    OtherAndUnknown =  2# (2,1,NA)value to use for cases where fire type is "OTHER" or "UNKNOWN",1 ="BURN",2="BUSHFIRE",NA = Fire excluded from analysis default is 2 ("BUSHFIRE")
                    
){
  myDF<-st_read(flattenedFH)[,c( "SEASON", "FIRETYPE")]
  names(myDF)[1]<- "SEASON"
  
  min.SEASON<-sort(unique(myDF$SEASON))[2]# second season in fire history
  
  #timespan ( range of consecutive years) for which Fire History sequences are calculated
  if((is.null(start.SEASON))|(start.SEASON<min.SEASON)){
    start.SEASON=min.SEASON
  } 
  
  if(is.null(end.SEASON)){
    max.SEASON<-max(myDF$SEASON)
  } else {
    max.SEASON = end.SEASON
  }
  
  TimeSpan<-start.SEASON:max.SEASON
  
  
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="BURN"]<-1
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="BUSHFIRE"]<-2
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="OTHER"]<-OtherAndUnknown
  myDF$FIRETYPE_NO[myDF$FIRETYPE=="UNKNOWN"]<-OtherAndUnknown
  
  myDF<-add_xystring(myDF)
  
  myDF<-unique(myDF[,c("geometry",
                       "XYString",
                       "SEASON",
                       "FIRETYPE_NO")])
  
  myDF<-myDF[with(myDF,order(XYString, SEASON, FIRETYPE_NO)),] 
  
  myDF$Sequence<-1 
  myDF$Sequence<-unlist(lapply(split(myDF$Sequence, myDF$XYString), cumsum)) 
  
  
  print("Making wide format firetype")
  TypeDF<-spread(myDF[,c("XYString","Sequence","FIRETYPE_NO")],Sequence,FIRETYPE_NO)
  ReLabelColNo<-2:(ncol(TypeDF)-1)
  names(TypeDF)[ReLabelColNo]<-paste("FireType",sprintf("%02d",as.numeric(names(TypeDF)[ReLabelColNo])),sep="")
  TypeDF[is.na(TypeDF)]<-0
  TypeDF<-TypeDF[,-1]
  TypeDFNoGeom<-TypeDF
  st_geometry(TypeDFNoGeom)<-NULL
  FT<-as.matrix(TypeDFNoGeom)
  
  
  print("Making wide format season")
  YearDF<-spread(myDF[,c("XYString","Sequence",'SEASON')],Sequence,SEASON)
  SEASNames<-paste("SEAS",sprintf("%02d",as.numeric(names(YearDF)[ReLabelColNo])),sep="")
  names(YearDF)[ReLabelColNo]<-SEASNames
  YearDF[is.na(YearDF)]<-0
  YearDF<-YearDF[,-1]
  YearDFNoGeom<-YearDF
  st_geometry(YearDFNoGeom)<-NULL
  M<-as.matrix(YearDFNoGeom)
  
  
  
  print ("calculating inter-fire intervals")
  Cols<-ncol(M)
  M[M==0]<-NA
  Interval<-M[,2:Cols]-M[,1:Cols-1]
  IntNames<-paste("INT",sprintf("%02d",1:(Cols-1)),sep="")
  colnames(Interval)<-IntNames
  
  OutDF<-cbind(YearDF,Interval)
  
  
  print("spatial join of wide format polygon datasets")
  OutDF<-st_join(OutDF,TypeDF,join=st_equals)
  #OutDF<-st_parallel(c(OutDF,TypeDF),st_join,join=st_equals,n_cores=4)
  OutTab<-OutDF
  st_geometry(OutTab)<-NULL
  
  
  Year<-as.matrix(OutTab[,SEASNames])
  
  
  print("making multicore cluster for parallel processing of TSF")
  #Ncores<-8#detectCores()-2
  cl<-makeCluster(Ncores)
  registerDoParallel(cl, cores=Ncores)
  
  print("made cluster, now calcuating time since fire for each input year")
  
  LBY<-foreach(y=iter(TimeSpan),.combine = cbind,.packages ="Rfast" )%dopar%
    try({
      X<-Year
      #X[X>=y]<-NA changed this line to deal with year shift issue.
      X[X>y]<-NA
      X[X==0]<-NA
      LBY<-rowMaxs(X,value=T)
      #LBY<-apply(X,1,max,na.rm=TRUE)
      LBY[is.infinite(LBY)]<-NA
      LBY
      
      
    })
  stopCluster(cl)
  tYSF<-TimeSpan-t(LBY)
  YSF<-t(tYSF)
  YSFNames<-paste0("YSF",TimeSpan)
  LBYNames<-paste0("LBY",TimeSpan)
  LFTNames<-paste0("LFT",TimeSpan)
  colnames(YSF)<-YSFNames
  colnames(LBY)<-LBYNames
  
  
  print("calculating lookup matrix for getting last firetype by year")
  LUM<-matrix(NA,nrow(M),max.SEASON)
  for (i in 1:nrow(M)){
    R<-i
    C<-as.numeric(na.omit(M[i,]))
    V<-(FT[i,(1:length(C))])
    LUM[R,C]<-V
  }
  LFT<-matrix(NA,nrow(M),length(TimeSpan))
  for(i in 1:nrow(M)){
    LFT[i,]<-LUM[i,LBY[i,]]}
  colnames(LFT)<-LFTNames
  
  OutDF<-cbind(OutDF,YSF)
  OutDF<-cbind(OutDF,LBY)
  OutDF<-cbind(OutDF,LFT)
  OutDF$ID<-as.integer(rownames(OutDF))
  print("completed making FH object")
  
  
  # FH_outpath<-outputFH
  # st_write(OutDF,outputFH)
  # print("written FH to disk")
  # 
  # 
  # 
  # 
  # # Using system command and gdal_rasterize since it is quicker than the R raster:rasterize 
  # 
  # outRaster <-FH_ID.tif
  # 
  # GDALRasterize(inVector=outputFH,outRaster=outRaster,ID<-"ID",TemplateRaster=CropDetails$Raster,datatype="INT4S")
  
  
  
  
  rm(ReLabelColNo,LUM,TypeDF,YearDF,TypeDFNoGeom,YearDFNoGeom)
  gc()
  results<-list(OutDF=OutDF,TimeSpan=TimeSpan,  YSFNames=YSFNames,LBYNames=LBYNames, LFTNames=LFTNames)
  return(results)
}

###makeYSF_LFT_YEAR_RASTERS#########################################################################################
#function to export rasters of TSF and YSF for all years # optional rest of script does not depend on these being made.

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
#make matrix of cell wise values of YSF and LFT 
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



###makeLU_List#############################################################################################
# function creates a list of Lookup array for each taxon(VBA_CODE) for  YSF x EFGNO x FireTypeNo
makeLU_List<- function(myHDMSpp_NO=HDMSpp_NO,myAbundDataLong=ExpertDataLong){
  
  
  
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

###calcDeltaAbund#######################################################################################

# Calculates baseline RA based on input baseline years and deviation from baseline for each future year, output written to CSV files.
#"SppSummChangeRelativetoBaseline.csv" and "SpYearSummSpreadbyYear.csv"

calcDeltaAbund<- function(SpYearSumm,
                          TimeSpan =rv$FHanalysis$TimeSpan,
                          myBaseline,
                          ResultsDir,
                          HDMSpp_NO,
                          TaxonList)
{
  SpYearSummSpreadbyYear<<-spread(SpYearSumm,key = SEASON,value = SUM_RAx100)
  write.csv(SpYearSummSpreadbyYear,file.path(ResultsDir,"SpYearSummSpreadbyYear.csv"),row.names=F)
  write.csv(TaxonList[TaxonList$TAXON_ID%in%HDMSpp_NO,],
            file.path(ResultsDir,"SppConsideredInAnalysis.csv"),
            row.names=F)
  
  #to get % of baseline need to define which columns provide the baseline ( one or mean of several using apply (mean)) then divide remaining values by this column.
  if (length(myBaseline==1)){
    Baseline<-SpYearSummSpreadbyYear[,as.character(myBaseline)]
  } else {
    Baseline<-apply(SpYearSummSpreadbyYear[,as.character(myBaseline)],1,mean)
  }
  SpYearSummSpreadbyYear$Baseline<-Baseline
  # gets the integer value for current year used in next line so that changes to baseline
  #are only displayed for future years or if no future yeas then years since baseline.
  ThisYear<-as.integer(format(Sys.Date(), "%Y"))
  SinceYear<-ifelse(sum(TimeSpan>ThisYear)>0,
                    ThisYear,
                    max(myBaseline))
  
  Deltas<-as.matrix(SpYearSummSpreadbyYear[,as.character(TimeSpan[TimeSpan>SinceYear])]/Baseline)
  names(Deltas)<-paste(names(Deltas),"prop baseline")
  NoLessthanThreshhold<-rowSums(Deltas<=SpYearSummSpreadbyYear$CombThreshold)
  LastLessThanThreshold<-Deltas[,ncol(Deltas)]<=SpYearSummSpreadbyYear$CombThreshold
  ChangeRelativeToBaseline<-cbind(SpYearSummSpreadbyYear[,c("TAXON_ID",
                                                            "ShortName",
                                                            "COMMON_NAME",
                                                            "NAME",
                                                            "DIVNAME",
                                                            "EPBC_ACT_STATUS",
                                                            "VIC_ADVISORY_STATUS",
                                                            "CombThreshold")],
                                  Deltas,
                                  NoLessthanThreshhold,
                                  LastLessThanThreshold)
  write.csv(ChangeRelativeToBaseline,
            file.path(ResultsDir,"SppSummChangeRelativetoBaseline.csv"),
            row.names=F)
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

###calcSppEFGLMU--------------------------------------------------------------------------------------------------------

#Calculate the species in each EFG in given area for GSO calcuations. works by using the indices of the standard dimesions 
#raster that are in the supplied shapefile region boundary ( via function maekCropDetails)
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

#get values of inputrasters depending on cellSize
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
#st_parallel--------------------------------------------------------
#possible parallel version of st_ functions  not yet expolored from
# https://www.spatialanalytics.co.nz/post/2017/09/11/a-parallel-function-for-spatial-analysis-in-r/ 
#and https://www.spatialanalytics.co.nz/post/2018/04/01/fixing-st-par/
#for linux only

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

cellsToHectares<-function(RasterMetres=RasterRes){
  (RasterMetres/100)^2
}
############END OF FUNCTIONS########################################
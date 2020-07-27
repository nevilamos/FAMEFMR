# calculate date sequence BBTFI (accomodating  Hi and Lo fire intesity of first burn to determine TFI) ----------------
#returns list containing  [[1]] the date sequence matrix for each cell of the raster
#                         [[2]] the EFG TFI Lookup for each cell of the raster
#                         [[3]] the raster resolution used.
# throws an error if raster resolution does not match selected FH_ID file

calc_BBTFI<-function(FHanalysis,#the slected FHanalysis object ( either through running analysis previously, or loading the rdata object.)
                     cropRasters,
                     TFI_LUT_DF = TFI_LUT#the dataframe read from a csv that gives the lookup table from EVD to MinTFI_LO, MinTFI_HI and MaxTFI
) {
  Hectares<-(as.numeric(FHanalysis$RasterRes)/100)^2
  
  r<-FHanalysis$FH_IDr
  
  FH_ID<-values(r)
  r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
            
            
            
            
            # FH_ID<-as.data.frame(values(r))
            # names(FH_ID)<-"ID"
            # DELWP<-cropRasters$DELWP
            # PLM<-cropRasters$PLM
            # EFG<-cropRasters$EFG
            # FIRE_REG<-cropRasters$RGN
            # FIREFMZ<-cropRasters$FIREFMZ
            # PLM<-cropRasters$PLM
            # EFG_DF<-as.data.frame(EFG)
            
            # TFI<-left_join(EFG_DF,TFI_LUT_DF)
            names(FH_ID)<-"ID"
            
            PLM<-cropRasters$PLM
            EFG<-cropRasters$EFG
            FIRE_REG<-as.integer(cropRasters$RGN)
            FIREFMZ<-as.integer(cropRasters$FIREFMZ)
            DELWP<-cropRasters$DELWP
            
            AllCombs<-as.data.table(cbind(FH_ID,EFG,FIRE_REG,FIREFMZ,PLM,DELWP))
            
            rm(PLM,EFG,FIRE_REG,FIREFMZ,FH_ID,DELWP)
            rm(cropRasters)
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
            TYPE_HI<-TYPE_HI[U_AllCombs_TFI$FH_ID,]
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
            #count of the number of times bbtfi events in each comination of EFG with fire history and grouping polygons
            TimesBBTFI<-rowSums(!is.na(BBTFI_COMB))
            TimesBBTFI[TimesBBTFI==0]<-NA
            FirstBBTFI<-apply(BBTFI_COMB,1,min,na.rm=T )
            FirstBBTFI[is.infinite(FirstBBTFI)]<-NA
            BBTFI_COMB<-cbind(U_AllCombs,BBTFI_COMB,TimesBBTFI,FirstBBTFI)
            
            
            #pivot tables to get one line per sequence eventx combination this may not be needed
            #BBTFI_COMB_LONG<-as.data.table(cbind(U_AllCombs,BBTFI_COMB))%>%pivot_longer(cols=all_of(SEASFields[-1]),names_to="SEAS_SEQ",values_to="SEASON")#there are only dates in season where it was BBTFI in that season


            
           
                                           
            
                                                                                       
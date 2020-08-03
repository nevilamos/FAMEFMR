calcU_All_Combs<- function(FHAnalysis,cropRasters){

TimeRange<-as.integer(FHanalysis$TimeSpan)
TimeNames<-as.character(FHanalysis$TimeSpan)
LTR<-length(TimeRange)
r<-FHanalysis$FH_IDr
FH_ID<-values(r)
r<-raster(nrows=nrow(r), ncols=ncol(r),ext=extent(r),crs=crs(r),vals=NULL)
gc()

#names(FH_ID)<-"ID"

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
U_AllCombs_TFI<-Join_Names(U_AllCombs_TFI)
U_AllCombs_TFI$Hectares=U_AllCombs_TFI$nPixel*cellsToHectares()
return(list("U_AllCombs_TFI"=U_AllCombs_TFI,"Index_AllCombs"=Index_AllCombs))
}

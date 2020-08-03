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

 #############################################################
makeGS_Sum<-function(writeGSRasters,
                     ResultsDir,
                     U_AllCombs_TFI =myAllCombs$U_AllCombs_TFI,
                     Index_AllCombs=myAllCombs$Index_AllCombs,
                     FHanalysis,
                     writeYears=NULL) {
  TimeNames<-as.character(FHanalysis$TimeSpan)
  GS_LU<-makeGS_LU()
  #read the corresponding shapefile for the FireHat analysis chosen
  #convert to a data.frame
  OutTab<-FHanalysis$OutDF
  st_geometry(OutTab)<-NULL
  #ID<-OutTab$ID
  
  YSF_Fields<-names(OutTab)[grep("^YSF",names(OutTab))]
  YSFplus1<-OutTab[,YSF_Fields]+1
  #partial inflation using U_AllCombs_TFI$FH_ID
  YSFplus1<-YSFplus1[U_AllCombs_TFI$FH_ID,]
  YSFplus1<-as.matrix(YSFplus1)
  #inflate GS_LU on EFG axis using U_AllCombs_TFI$$EFG
  GS_LU<-t(GS_LU)[U_AllCombs_TFI$EFG,]

  indmat=cbind(rep(1:nrow(YSFplus1),ncol(YSFplus1)),as.vector(YSFplus1))
  

  GSrows=matrix(GS_LU[indmat],nrow=nrow(YSFplus1))
  
  colnames(GSrows)<-TimeNames

  GS_Summary<-cbind(U_AllCombs_TFI,GSrows)
  
  
  GS_Summary_Long<-GS_Summary%>%
    select(-c( MIN_LO_TFI,MIN_HI_TFI,MAX_TFI,Index,FH_ID,FIRE_REG,FIREFMZ,DELWP ))%>%
    pivot_longer(all_of(TimeNames),names_to="SEASON",values_to="GS")%>%
    group_by(EFG,EFG_NAME,PLM ,FIRE_FMZ_NAME,FIRE_FMZ_SHORT_NAME,FIRE_REGION_NAME,DELWP_REGION,SEASON,GS)%>%
    summarise(Pixels=sum(nPixel),Hectares=sum(Hectares))
  GS_Summary_Long<-left_join(GS_Summary_Long,GS_LUT)
  
 
  setnames(GS_Summary,old=TimeNames,new=GS_Names<-paste0("GS_",TimeNames))
 return(list("GS_Summary_wide"=GS_Summary, "GS_Summary_Long"= GS_Summary_Long))
}
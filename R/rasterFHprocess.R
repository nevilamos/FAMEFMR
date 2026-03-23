
#' Raster FireHisotry Process
#'
#' @param fireHistory
#' @param template_raster
#' @param OtherAndUnknown
#'
#' @return
#' @export
#'
#' @examples
raster_fh_process<-function(fireHistory=mySF, template_raster=template_r,OtherAndUnknown=2){

  mySF<-fireHistory
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BURN"] <- 1
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "BUSHFIRE"] <- 2
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "OTHER"] <- OtherAndUnknown
  mySF$FIRETYPE_NO[mySF$FIRETYPE == "UNKNOWN"] <- OtherAndUnknown

  mySF$FireDetail<-as.integer(paste0(mySF$SEASON,mySF$FIRETYPE_NO))

  LayerNames=sort(unique(mySF$FireDetail))
  L<-length(LayerNames)

  SEASONr<-terra::rast(template_raster,nlyrs=L,names=LayerNames)
  FTr<-terra::rast(template_raster,nlyrs=L,names=LayerNames)
  #have to use 0 rather than NA otherwise rows that are all NA excluded from terra::as.data.frame()
  for(i in LayerNames){
    SEASONr[[as.character(i)]]<-terra::rasterize(mySF %>% dplyr::filter(FireDetail==i),y = template_raster,background=0,field = "SEASON")
    FTr[[as.character(i)]]<-rasterize(mySF %>% filter(FireDetail==i),y = template_raster,background=0,field = "FIRETYPE_NO")
  }

  SEASONm <-  values(SEASONr,mat=T)
  shift_zero_in_place(SEASONm)
  SEASONm<-drop_zero_cols(SEASONm)
  uSEASONm<-unique(SEASONm)
  SEASNames<-paste0("SEAS",sprintf("%02d",1:ncol(uSEASONm)))
  colnames(uSEASONm)<-SEASNames
  uSEASONdf<-as_tibble(uSEASONm) %>% mutate(SEAS_ID =row_number(),.before=1)




  FTm <-  values(FTr,mat=T)
  FTm<-shift_zero(FTm)# for some reason this less memory efficient function works here but shift zero in place does not
  FTm<-drop_zero_cols(FTm)
  uFTm<-unique(FTm)
  FTNames<-paste0("FireType",sprintf("%02d",1:ncol(uFTm)))
  colnames(uFTm)<-FTNames
  uFTdf<-as_tibble(uFTm) %>% mutate(FT_ID =row_number(),.before=1)



  r<-rast(template_raster,nlyrs=ncol(SEASONm))
  values(r)<-SEASONm
  names(r)<-SEASNames

  rSeasOut<-subst(r,unique(uSEASONm),1:nrow(uSEASONm))
  names(rSeasOut)<-"SEAS_ID"


  values(r)<-FTm
  names(r)<-FTNames

  rFTOut<-subst(r,uFTm,1:nrow(uFTm))
  names(rFTOut)<-"FT_ID"

  rFT_SEAS<-c(rSeasOut,rFTOut)
  FT_SEASm<-values(rFT_SEAS,matrix=T)
  uFT_SEASm<-unique(FT_SEASm)
  uFT_SEASdf<-as_tibble(uFT_SEASm) %>% mutate(FH_ID =row_number(),.before=1) %>% left_join(uSEASONdf) %>% left_join(uFTdf)
  FH_IDr<-subst(rFT_SEAS,uFT_SEASdf %>% select(SEAS_ID,FT_ID),uFT_SEASdf %>% select(FH_ID) %>% pull())

  names(FH_IDr)<-"FH_ID"
  return(list(FH_IDr,uFT_SEASdf))

}

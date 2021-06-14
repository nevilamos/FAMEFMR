#' Title Summary of JFMP areas burned and not burned
#'
#' @param JFMP_DF JFMP table with columns for FMZ_ID,PU,and JFMPStatus plus all the scores calculated in jfmp1()
#'
#' @return JFMP_Summary data frame CSV reporting table 2, with one row for each District and columns for:
#'– Hectares allocated to burns in auto-JFMP in each zone
#'– Total hectares allocated to burn
#'– Score for each metric (x4) if JFMP implemented
#'– Score for each metrics (x4) if JFMP not implemented
#'- Hectares allocated in each FMZ_ID
#' @export
#'
#' @examples
jfmpSummary<- function(JFMP_DF="JFMPSummary with burned or unburned"){

JFMP_DF = AutomJFMP_DF
  JFMP_Summ<-JFMP_DF%>%
    ungroup()%>%
    mutate(FaunaRAImp=if_else(AutoJFMP_State=="BURN",WtSumRA_Burn,WtSumRA_NoBurn))%>%
    mutate(NeverBBTFI_Imp=if_else(AutoJFMP_State=="BURN",PuAreaHa- (ifelse(is.na(Burn_BBTFI),0,Burn_BBTFI)) - (ifelse(is.na(NoBurn_BBTFI),0,NoBurn_BBTFI)),PuAreaHa- (ifelse(is.na(NoBurn_BBTFI),0,NoBurn_BBTFI))))%>%
    mutate(LP1_Imp =if_else(AutoJFMP_State=="BURN",LP1_Burn,LP1_NoBurn))%>%
    mutate(LP2_Imp =if_else(AutoJFMP_State=="BURN",LP2_Burn,LP2_NoBurn))%>%
    group_by(DISTRICT_N)%>%
    summarize(AreaHa =sum(PuAreaHa,na.rm=T),
              SumFaunaRAImp=sum(FaunaRAImp,na.rm=T),
              SumFaunaRANoJFMP=sum(WtSumRA_NoBurn,na.rm=T),
              SumNeverBBTFI_Imp=sum(NeverBBTFI_Imp,na.rm=T),
              SumNeverBBTFI_NoJFMP=sum(NoBurn_BBTFI,na.rm=T),
              SumLP1_Imp=sum(LP1_Imp,na.rm=T),
              SumLP1_NoJFMP=sum(LP1_NoBurn,na.rm=T),
              SumLP2_Imp=sum(LP2_Imp,na.rm=T),
              SumLP2_NoJFMP=sum(LP2_NoBurn,na.rm=T))
  
  AreasBurnedUnburned<-JFMP_DF%>%
    ungroup()%>%
    group_by(DISTRICT_N,AutoJFMP_State)%>%
    summarise(AreaHa=sum(PuAreaHa,na.rm=T))%>%
    spread(AutoJFMP_State,AreaHa)%>%
    rename(Burned_ha='BURN',NoBurn_Ha='NO BURN')
  
  AreasBurnedbyFMZ_CODE<-JFMP_DF%>%
    ungroup()%>%
    filter(AutoJFMP_State=='BURN')%>%
    group_by(DISTRICT_N,X_ZONETYPE)%>%
    summarise(AreaHa=sum(PuAreaHa,na.rm=T))%>%
    spread(X_ZONETYPE,AreaHa)
  names(AreasBurnedbyFMZ_CODE)[-1]<-paste("Burned_ha ",names(AreasBurnedbyFMZ_CODE)[-1])
    
  JFMP_Summ<-left_join(AreasBurnedUnburned,JFMP_Summ)
  JFMP_Summ<-left_join(JFMP_Summ,AreasBurnedbyFMZ_CODE)
  return(JFMP_Summ)
  
}

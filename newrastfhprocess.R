library(tidyterra)
library(terra)
library(data.table)
library(Rcpp)



myFH<-vect("rawFH/rawFH_FIRE_HISTORY_LF_DISTRICT_1755.gpkg")
myFH$FIRETYPE_NO[myFH$FIRETYPE == "BURN"] <- 1
myFH$FIRETYPE_NO[myFH$FIRETYPE == "BUSHFIRE"] <- 2
myFH$FIRETYPE_NO[myFH$FIRETYPE == "OTHER"] <- 2
myFH$FIRETYPE_NO[myFH$FIRETYPE == "UNKNOWN"] <- 2

myFH<-myFH%>% mutate(SEASON_TYPE =(SEASON*100)+FIRETYPE_NO)

SeasonTypes<-sort(unique(myFH$SEASON_TYPE))
r<-terra::rast("InputGeneralRasters/EFG_NUM_225.tif")
crs(r)<-crs(myFH)

s<-NULL
for (i in SeasonTypes){
x<-rasterize(myFH %>% filter(SEASON_TYPE==i),r)
names(x)<-as.character(i)
s<-c(s,x)

}
rm(x)
rm(i)


s<-rast(s)
v<-values(s)

out <- extract_season_firetypes_with_group_id(v, colnames(v))
rm(s)
rm(v)
gc()
# Convert to data.tables
out$SEASONS <- as.data.table(out$SEASONS)
out$FireTypeNos <- as.data.table(out$FireTypeNos)


# Insert group_id as the first column
out$SEASONS[, group_id := out$group_id]
setcolorder(out$SEASONS, c("group_id", setdiff(names(out$SEASONS), "group_id")))

FireTypeNos[, group_id := out$group_id]
setcolorder(FireTypeNos, c("group_id", setdiff(names(FireTypeNos), "group_id")))



uSEASONS<-leftPackMissingOMP(as.matrix(unique(out$SEASONS)))
uFireTypeNos<-leftPackMissingOMP(as.matrix(unique(out$FireTypeNos)))




# SEASONS<-as.data.table(t(t(v)*as.integer(substr(colnames(v),1,4))))
# FireTypeNos<-as.data.table(t(t(v)*as.integer(substr(colnames(v),6,6))))
# rm(s)
# rm(v)
# gc()
# SEASONS[, group_id := .GRP, by = names(SEASONS)]
# uSEASONS<-unique(SEASONS)
# FireTypeNos[, group_id := .GRP, by = names(FireTypeNos)]
# uFireTypeNos<-unique(FireTypeNos)
#
#
#
# removeNaN<-function(x){
#   x<-x[is.finite(x)]
#   L=length(x)
#   x<-as_tibble(t(x))
#   if(L>0)
#   names(x)<-as.character(1:L)
#   return(x)
#   }
#
#
#
#
#
#
# # Number of cores to use
# no_cores <- parallel::detectCores() - 1
#
# # Create the cluster
# cl <- makeCluster(no_cores)
# registerDoSNOW(cl)
#
# # Run parallel foreach loop over rows
# result_list <- foreach::foreach(i = 1:nrow(SEASONS), .combine = 'bind_rows') %dopar% {
#   removeNaN(SEASONS[i, ])
# }
#
# # Stop cluster
# stopCluster(cl)
#
# # Bind results together
# result_df <- bind_rows(result_list)
# uSEASONS1<-bind_rows(apply(SEASONS,MARGIN = 1,FUN = removeNaN))
# lapply(uSEASONS1,bind_rows)
# unlist(lapply(uSEASONS1,length))

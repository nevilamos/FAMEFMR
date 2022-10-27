# #r prperpocessing working with to y data but error on real dataset Error: [union] geom types do not match even when weeding out multipolygons
#
# library(terra)
# library(sf)
# sff<-vect("D:/FAMEshiny/FAMEPreProcessing/commondata/SimpleFireFuture.shp")
# sfh<-vect("D:/FAMEshiny/FAMEPreProcessing/commondata/SimpleFireHistory.shp")
# x<-vect("D:/FAMEshiny/FAMEPreProcessing/commondata/FIRE_HISTORY_20210310.shp")
# x<-st_as_sf(x)
# x<-x[st_geometry_type(x)== "POLYGON",]
# x<-st_buffer(x,0)
# unique(st_geometry_type(x))
# x<-vect(x)
#
#
# #the below x works
# x<-rbind(sff,sfh)
#
# x$x_row<-1:nrow(x)
# x_df<-as.data.frame(x)
#
# y<-terra::union(x)
#
# y<-st_as_sf(y)
# y_long<-  y %>%
#   tidyr::pivot_longer(-geometry) %>%
#   mutate(x_row = parse_number(name)) %>%
#   filter(value == 1) %>%
#   dplyr::left_join(x_df)

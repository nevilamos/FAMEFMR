library(FAMEFMR)
library(terra)
# nFires<-20 #number of fires to generate in dummy dataset
# inFH <- make_organic_polys(seed=4,n = nFires)
# inFH$SEASON<-sample(1939L:2025L,replace = T,size = nFires)
# inFH$FIRETYPE<-sample(c("BUSHFIRE","BURN"),replace = T,size = nFires)
# inFH<-inFH[c("SEASON","FIRETYPE")]
infile<-"../FAMEshiny/rawFH/demoFH2025_vg94.gpkg"#"../FAMEshiny/rawFH/rawFH_fh_rep_01_Central_highlands.gpkg"


# writeVector(inFH,infile,overwrite=TRUE)
v<-vect(infile)
template_r<-rast(ext=ext(v),res=225,crs=crs(v))
rm(v)
# r_template = template_r
# vector_file = infile
# fields = c("SEASON", "FIRETYPE_NO")
# combField_multiplier = 10L
# background = 0
# overwrite = TRUE
# datatype = "INT2U"
# compress = "LZW"
# memfrac = .8
# progress = TRUE
# OtherAndUnknown = 2
# start.SEASON = 1980
# end.SEASON =NA
# max_interval=0





fhAnalysis<-fhProcess_raster_first()




#loading dev versions for testing
# source("C:/Data/FAMEFMR/R/make_fh_raster_stack_dev.R", echo=TRUE)
# source("C:/Data/FAMEFMR/R/process_blocks_to_bin_dev.R", echo=TRUE)
# source("C:/Data/FAMEFMR/R/write_raster_from_bin_dev.R", echo=TRUE)
# source("C:/Data/FAMEFMR/R/add_fire_lft_lby_ysf_dev.R", echo=TRUE)
#fhProcess_raster_first()
s<-make_fh_raster_stack()

#r<-rast(dir(out_dir,full.names = T))

p1<-process_blocks_to_bin(s)

r2<-write_raster_from_bin(r_template = r_template,bin_file = p1$bin,index_file = p1$index,max_ncol = p1$max_ncol)

um<-as.matrix(unique(r2))
FH_ID<-1:nrow(um)
FH_IDr<-subst(x = r2,um,FH_ID)

ums<-split_integer_first_last(as.integer(um),5,4,1 )


FTm<-SEASm<-um
SEASm[]<-ums$first
colnames(SEASm)<-paste0("SEAS",sprintf("%02d",1:ncol(um)))
FTm[]<-ums$last
colnames(FTm)<-paste0("FireType",sprintf("%02d",1:ncol(um)))

OutDF<-cbind(FH_ID,SEASm,FTm)
FH_IDr
fhAnalysis<-add_fire_lft_lby_ysf (OutDF = OutDF,max_interval = max_interval,start.SEASON = start.SEASON,end.SEASON = end.SEASON,v=vect(infile))

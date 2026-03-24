infile <- "../FAMEshiny/rawFH/demoFH2025_vg94.gpkg"
v<-vect(infile)
template_r<-rast(ext=ext(v),res=2,crs=crs(v))
r_template = template_r
vector_file = infile
fields = c("SEASON", "FIRETYPE_NO")
combField_multiplier = 675L
background = 0
overwrite = TRUE
datatype = "INT2U"
compress = "LZW"
memfrac = .8
progress = TRUE
gc_every = 0
OtherAndUnknown = 2
start.SEASON = 1980
end.SEASON =NA
out_dir<-file.path(tempdir(),"FIRE_RASTER_STACK")


make_fh_raster_stack()

r<-rast(dir(out_dir,full.names = T))

p1<-process_blocks_to_bin(r)

write_raster_from_bin(r_template = r,bin_file = p1$bin,index_file = p1$index,max_ncol = p1$max_ncol,out_raster_file = "testout.tif")
r2<-rast("testout.tif")
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



SEASNames <- colnames(OutDF)[grep(pattern = "SEAS", colnames(OutDF))]
FTNames <- colnames(OutDF)[grep(pattern = "FireType", colnames(OutDF))]
SEAS_Matrix <- as.matrix(OutDF[, SEASNames])
FT_matrix <- as.matrix(OutDF[, FTNames])
Cols <- ncol(SEAS_Matrix)
SEAS_Matrix[SEAS_Matrix == 0] <- NA
Interval <- SEAS_Matrix[, 2:Cols] - SEAS_Matrix[, 1:Cols -
                                                  1]
IntNames <- paste("INT", sprintf("%02d", 1:(Cols - 1)), sep = "")
colnames(Interval) <- IntNames
OutDF <- cbind(OutDF, Interval)
if (max_interval > 0) {
  FT_matrix <-
    fireTypeLowToHigh(
      max_interval = as.integer(max_interval),
      Interval_Matrix = Interval,
      Firetype_Matrix = FT_matrix
    )
  OutDF[, FTNames] <- FT_matrix
}else if (max_interval < 0) {
  stop("max interfal cannot be less than 0")
}else {}



min.SEASON <- sort(unique(v$SEASON))[2]
if (is.na(start.SEASON)) {
  start.SEASON = min.SEASON
}else {
  if (start.SEASON < min.SEASON) {
    start.SEASON = min.SEASON
  }
  else {
    start.SEASON = start.SEASON
  }
}
if (is.na(end.SEASON)) {
  max.SEASON <- max(v$SEASON)
}else {
  #max.season must be greater than or equal to maximum seaon in FH ( including
  #JFMP0+2 case)
  max.SEASON = max(end.SEASON,v$SEASON)
}
TimeSpan <- start.SEASON:max.SEASON








LTR <- length(TimeSpan)
SEAS_Matrix[is.na(SEAS_Matrix)] <- 0
LBY <- matrix(NA, nrow(SEAS_Matrix), LTR)
for (i in 1:LTR) {
  try({
    y = TimeSpan[i]
    LBY[, i] <- LBY_f(M = SEAS_Matrix, y)
  })
}
tYSF <- TimeSpan - t(LBY)
YSF <- t(tYSF)
YSFNames <- paste0("YSF", TimeSpan)
LBYNames <- paste0("LBY", TimeSpan)
LFTNames <- paste0("LFT", TimeSpan)
colnames(YSF) <- YSFNames
colnames(LBY) <- LBYNames
print("calculating lookup matrix for getting last FireType by SEASON")
SEAS_Matrix[SEAS_Matrix == 0] <- NA
LUM <- matrix(NA, nrow(SEAS_Matrix), max.SEASON)
for (i in 1:nrow(SEAS_Matrix)) {
  R <- i
  C <- as.numeric(stats::na.omit(SEAS_Matrix[i,]))
  V <- (FT_matrix[i, (1:length(C))])
  LUM[R, C] <- V
}

print("calculating last fire type")
LFT <- matrix(NA, nrow(SEAS_Matrix), LTR)
for (i in 1:nrow(SEAS_Matrix)) {
  LFT[i,] <- LUM[i, LBY[i,]]
}
colnames(LFT) <- LFTNames
OutDF <- cbind(OutDF, YSF)
OutDF <- cbind(OutDF, LBY)
OutDF <- cbind(OutDF, LFT)

############################################################################
# This function  compiles vectors of raster values  stored in ./HDMS/75m or
# ./HDMS/225m (75 m or 225 m versions), to filematrix files on disk - these
# are large, but faster to read the data from than reading directly from
# raster or an .rdata file 

# It takes a fair while to run - but only needs to be  run if the rasters are
# updated, or new rasters added.  The rasters are  values are stored as a single
# column vector per raster with the character from of the VBA TAXON_ID as the
# column name - this is later used to retrieve the values in species
# calculations currently these rasters contain binary values - and act as a
# species  mask on the abundace values. The rasters could be updated to
# continuous or thresholded continuous values - if this is done change the
# commented line below from the one with type="logical" to type = "double"
# This will result in greater disk storage space requirments for the filematrix
# files, and also memory for cacluations

# written by nevil.amos@delwp.vic.gov.au
############################################################################

#individual new rasters values can be appended easily to the fm objects, or
#individual species rasters replaced by using the colname - this is faster than
#repeating the whole process ( which takes >1 hour for 577 75m rasters)

for (j in c(75, 225)){
  myTiffs <- dir(file.path(HDMDir, paste0(j, "m")), pattern = ".tif", recursive = TRUE, full.names = TRUE)
  myCols <- length(myTiffs)
  myRows <- ncell(raster(myTiffs[1]))
  
  fm <- fm.create(filenamebase = file.path("HDMS", (paste0("HDMVals", j))),
                  nrow = myRows,
                  ncol = myCols,
                  type = "logical")
  
  for (i in 1:myCols){
    fm[,i] <- values(raster(myTiffs[i]))
    print(i)
  }
  colnames(fm) <- as.character(get_Spp_No(myTiffs))
}

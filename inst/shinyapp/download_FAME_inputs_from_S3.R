
# library(aws.s3)
newDirPaths<-c("AdHocPolygons",
               "CustomCSV",
               "FH_Outputs",
               "HDMS",
               "InputGeneralRasters",
               "PUPolygons",
               "rawFH",
               "ReferenceShapefiles",
               "ReferenceTables",
               "results"
)

for(mypath in newDirPaths){
  
  if (file.exists(mypath)){
    
  } else {
    dir.create(mypath,recursive = T)
    
  }
  #   
  #   
  
}
datazip<-"FAME4_VG2020_inputfiles.zip"
if(!file.exists(datazip)){
  download.file("https://ecological-risk-analysis.s3.ap-southeast-2.amazonaws.com/FAME_FMR/FAME4_VG2020_inputfiles.zip",datazip)
  
  if (.Platform$OS.type == "unix"){
    system("unzip -o ./FAME4_VG2020_inputfiles.zip && cp -r FAME4_VG2020_inputfiles/. . ")
  }else{
    unzip("./FAME4_VG2020_inputfiles.zip",exdir = ".")
    system("xcopy .\\FAME4_inputfiles\ . /E/H/y")
    
  }
  unlink(file_path_sans_ext(datazip),recursive=TRUE)
}




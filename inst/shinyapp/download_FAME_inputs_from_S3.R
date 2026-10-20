
# library(aws.s3)
newDirPaths<-c(fame_data_path("AdHocPolygons"),
               fame_data_path("CustomCSV"),
               fame_data_path("FH_Outputs"),
               fame_data_path("HDMS"),
               fame_data_path("InputGeneralRasters"),
               "PUPolygons",
               fame_data_path("rawFH"),
               fame_data_path("ReferenceShapefiles"),
               fame_data_path("ReferenceTables"),
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




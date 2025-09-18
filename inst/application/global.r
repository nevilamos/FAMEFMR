rm(list = ls(all = TRUE))


options(stringsAsFactors = F)
source("installationCheck.R")
# get version of FAMEFMR in use and set app version

FAMEGUIVersion<-FAMEFMRVersion <- paste("  R", getRversion(), "FAMEFMR", packageVersion("FAMEFMR"))

#set option so that readr does not print column types
options(readr.show_col_types = FALSE)

# file that disables buttons while processes are running to prevent multi-clicks
# and lengthy operations running multiple times
# source("ButtonDisableHelpers.r")



# Set the maximum size of files for upload/ download

options(shiny.maxRequestSize = 5 * 1024^3)


# code to implement action on closing browser window
jscode <- "shinyjs.closeWindow = function() { window.close(); }"
####################################################################################

# MAKE RESULTS DIRECTORIES
# create a unique results directory for each run of scenarios
# using starting time as a string for Results directory name
# this is zipped for download of results.


StartTimeString <- format(Sys.time(), "%Y%m%d_%H%M")

WD <- getwd()

#cleans up any old empty results directories
if(.Platform$OS.type == "unix") {
  system("find ./results -type d -empty -delete")
  if (!dir.exists("./results")){
    dir.create("results")
  }

} else {
  #not currently working for windows
}

# Makes resultsDir
resultsDir <- file.path("./results", StartTimeString)

#
for (d in c(resultsDir)) {
  dir.create(d)
}
rm(d)
dir.create(file.path(resultsDir, "RA_Rasters"))
dir.create(file.path(resultsDir, "TFI_Rasters"))


# make lookup tables used in app
source("makeLUTS.R")

cellSizes <- c("225" = 225, "75" = 75)
#HDMVals225 <- "./HDMS/HDMVals225.rdata"
writeSpRasters <- FALSE # TRUE if rasters are to be output ( large number of files and considerable disk space)

# # Lookup table from EFG to TFI attributes ( csv version of CGDL lookup table)
# TFI_LUT <-
#   read_csv("./ReferenceTables/EFG_EVD_TFI.csv")[, c("EFG_NUM", "MIN_LO_TFI", "MIN_HI_TFI", "MAX_TFI", "EFG_NAME")]
# names(TFI_LUT)[1] <- "EFG"
#
#
# EFG_TSF_4GS <-
#   read_csv("./ReferenceTables/EFG_TSF_4GScorrectedAllEFGto400yrsV2.csv")[, c("EFG_NO", "GS4_NO", "YSF")]
#
# outputNames<-read.csv("./ReferenceTables/outputNames.csv",row.names =1)

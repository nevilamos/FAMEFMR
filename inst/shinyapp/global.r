#rm(list = ls(all = TRUE))




options(shiny.sanitize.errors = FALSE)
options(shiny.error = browser)
options(shiny.fullstacktrace = TRUE)

# --- CRASH DEBUGGER START ---
trace(shiny:::applyInputHandlers, at = 1, tracer = quote({
  cat("\n=== Incoming Input IDs and Types ===\n")
  for (n in names(inputs)) {
    cat(sprintf("Input ID: %s | Type: %s\n", n, attr(inputs[[n]], "type")))
  }
  cat("====================================\n\n")
}))
# --- CRASH DEBUGGER END ---

# library(reactlog)
# reactlog_enable()
## If a package is installed, it will be loaded. If any
## are not, the missing package(s) will be installed
## from CRAN and then loaded.

## First specify the packages of interest
packages = c("aws.s3",
             "dplyr",
             "ggplot2",
             "knitr",
             "nloptr",
             "plotly",
             "qs",
             "readr",
             "readxl",
             "Rfast",
             "rmarkdown",
             "rlang",
             "scales",
             "sf",
             "shiny",
             "shinycssloaders",
             "shinydashboard",
             "shinyFiles",
             "shinyjs",
             "shinyWidgets",
             "stringr",
             "terra",
             "tibble",
             "tidyr",
             "tools")

## Now load or install&load all packages from cran
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


options(stringsAsFactors = F)
#source("installationCheck.R")
# get version of FAMEFMR in use and set app version

#FAMEGUIVersion <- "Version 4.6 April 30 2026 VicGrid2020"

FAMEFMRVersion <- paste("  R", getRversion(), "FAMEFMR", packageVersion("FAMEFMR"))

#set option so that readr does not print column types
options(readr.show_col_types = FALSE)

# file that disables buttons while processes are running to prevent multi-clicks
# and lengthy operations running multiple times
#source("ButtonDisableHelpers.r")



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
resultsDir <- fame_data_path("results", StartTimeString)

#
for (d in c(resultsDir)) {
  dir.create(d)
}
rm(d)
dir.create(file.path(resultsDir, "RA_Rasters"))
dir.create(file.path(resultsDir, "TFI_Rasters"))


# make lookup tables used in app
## MAKE LOOK UP TABLES

# list of options for choice of region/ state or adhoc polygon for analysis
#ultimately should be possible to remove this and use  REG_LUT for the same purposes.
REG_NO <- c("USER DEFINED POLYGON" = 7,
            "WHOLE OF STATE" = 99,
            "BARWON SOUTH WEST" = 1,
            "GIPPSLAND" = 2 ,
            "GRAMPIANS" = 3,
            "HUME" = 4,
            "LODDON MALLEE" = 5,
            "PORT PHILLIP" = 6
)

REG_LUT <- data.frame(
  FIRE_REG = c(7,99, 1, 2, 3, 4, 5, 6),
  FIRE_REGION_NAME = c(
    "USER DEFINED POLYGON",
    "WHOLE OF STATE",
    "BARWON SOUTH WEST",
    "GIPPSLAND",
    "GRAMPIANS",
    "HUME",
    "LODDON MALLEE",
    "PORT PHILLIP"
  )
)

# Fire FMZ look up table
FIREFMZ_LUT <- data.frame(
  FIREFMZ = c(0, 1, 2, 3, 4, 5),
  FIRE_FMZ_NAME = c(
    "0 - Non FPA Land",
    "1 - Asset Protection Zone",
    "2 - Bushfire Moderation Zone",
    "3 - Landscape Management Zone",
    "4 - Planned Burn Exclusion Zone",
    "5 - Unknown. Contact Fire Management Officer"
  ),
  FIRE_FMZ_SHORT_NAME = c("Non FPA",
                          "APZ",
                          "BMZ",
                          "LMZ",
                          "PBEZ",
                          "UNK")
)

# DELWP region names look up table
DELWP_LUT <-  data.frame(
  DELWP = c(2, 4, 3, 5, 6, 7),
  #####-----LUT and then you could just slice that LUT.
  DELWP_REGION = c(
    "GIPPSLAND",
    "HUME",
    "PORT PHILLIP",
    "BARWON SOUTH WEST",
    "LODDON MALLEE",
    "GRAMPIANS"
  )
)

# TFI_STATUS values lookup table and export to csv
TFI_STATUS_LUT <- structure(
  list(
    TFI_VAL = c(-99, 0, 1, 5, 6),
    TFI_STATUS = c(
      "NONE",
      "WITHIN_TFI",
      "BELOW_MIN_TFI",
      "ABOVE_MAX_TFI",
      "ABOVE_MAX_BELOW_MIN_HI_TFI"
    )
  ),
  row.names = c(NA, -5L),
  class = c("tbl_df", "tbl", "data.frame")
)


# Firetype look up table
FIRETYPE_LUT <- data.frame(TYPE = c(1, 2),
                           FIRETYPE = c("BURN", "BUSHFIRE"))

# Growth Stage look up table
GS_LUT <- data.frame(
  "GS" = c(0, 1, 2, 3, 4),
  "GROWTH_STAGE" = c("Unknown", "Juvenile", "Adolescent", "Mature", "Old")
)

# EFG to TFI attributes look up table
# read csv version of CGDL lookup table
TFI_LUT <-
  readr::read_csv(fame_data_path("ReferenceTables", "EFG_EVD_TFI.csv"))[, c("EFG_NUM", "MIN_LO_TFI",
                                                                            "MIN_HI_TFI", "MAX_TFI", "EFG_NAME")]
names(TFI_LUT)[1] <- "EFG"


cellSizes <- c("225" = 225, "75" = 75)
HDMVals225 <- fame_data_path("HDMS", "HDMVals225.rdata")
writeSpRasters <- FALSE # TRUE if rasters are to be output ( large number of files and considerable disk space)

# Lookup table from EFG to TFI attributes ( csv version of CGDL lookup table)
TFI_LUT <-
  readr::read_csv(fame_data_path("ReferenceTables", "EFG_EVD_TFI.csv"))[, c("EFG_NUM", "MIN_LO_TFI", "MIN_HI_TFI", "MAX_TFI", "EFG_NAME")]
names(TFI_LUT)[1] <- "EFG"


EFG_TSF_4GS <-
  read_csv(fame_data_path("ReferenceTables", "EFG_TSF_4GScorrectedAllEFGto400yrsV2.csv"),
           col_types = cols(.default = col_integer()) # forces numeric columns to be read as integers, needed to ensure that output raster attribute table for growth stages is also as integer.
  )[, c("EFG_NO", "GS4_NO", "YSF")]

outputNames<-read.csv(fame_data_path("ReferenceTables", "outputNames.csv"),row.names =1)

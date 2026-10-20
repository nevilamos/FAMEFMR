
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

#check that FAMEFMR package is installed and get it from github if not.
if(require("FAMEFMR")==FALSE){
  devtools::install_github("nevilamos/FAMEFMR@v1.1")
}


#run once only to copy all the input files from AWS bucket-----------------------------
source("download_FAME_inputs_from_S3.R")
#for debugging------------------
options(warn = -1)

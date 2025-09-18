
## If a package is installed, it will be loaded. If any
## are not, the missing package(s) will be installed
## from CRAN and then loaded.

## First specify the packages of interest
packages = c("aws.s3",
             "dplyr",
             "ggplot2",
             "knitr",
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
             "tools",
             "foreach",
             "parallel",
             "doParallel")

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

#for debugging------------------
options(warn = -1)


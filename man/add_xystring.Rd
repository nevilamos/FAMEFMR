% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{add_xystring}
\alias{add_xystring}
\title{Adds concatenated String of X and Y coordinates of centroids of polygons to
Simple Features polygon object. This String acts as a key to identify
spatially identical polygons for use in tidyverse pivot functions.}
\usage{
add_xystring(myDF)
}
\arguments{
\item{myDF}{sf polygon object}
}
\value{
character vector of XYStrings
}
\description{
Adds concatenated String of X and Y coordinates of centroids of polygons to
Simple Features polygon object. This String acts as a key to identify
spatially identical polygons for use in tidyverse pivot functions.
}

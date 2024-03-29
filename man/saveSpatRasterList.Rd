% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveSpatRasterList.R
\name{saveSpatRasterList}
\alias{saveSpatRasterList}
\title{Save List containing SpatRaster objects ( package terra)}
\usage{
saveSpatRasterList(
  myList,
  filePath = "filepath to save list file to ending in qs"
)
}
\arguments{
\item{myList}{a list object that may contain SpatRasters}

\item{filePath}{file path to save the list object to ( as a package qs file with extention qs)}
}
\description{
Save List containing SpatRaster objects ( package terra)
}
\details{
Allows saving and reloading of a list from and R session that contains SpatRasters the function recursively searches for SpatRaster objects and coverts them to "PackedSpatRaster" using terra::wrap()
A work around to allow saving of sessions of R containing SpatRasters that can be reloaded ( saving the SpatRaster itself results in a null pointer error)
}

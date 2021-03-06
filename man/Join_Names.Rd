% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{Join_Names}
\alias{Join_Names}
\title{Joins one or more lookup tables to table containing ID values Function joins
Lookup tables (LUTS) to dataframe containing ID_NO: Name combinations}
\usage{
Join_Names(myDF, LUTS = c("TFI_LUT", "FIREFMZ_LUT", "REG_LUT", "DELWP_LUT"))
}
\arguments{
\item{myDF}{dataframe or similar containing indices for the LUTS listed, to
which the LUTS will be dplyr::left_joined}

\item{LUTS}{vector of names of LUTS in memory defaults
=c("TFI_LUT","FIREFMZ_LUT","REG_LUT","DELWP_LUT")}
}
\value{
a data.frame with the LUTS joined to it
}
\description{
Joins one or more lookup tables to table containing ID values Function joins
Lookup tables (LUTS) to dataframe containing ID_NO: Name combinations
}

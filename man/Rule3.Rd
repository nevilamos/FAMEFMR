% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSOFunc_Rule3.R
\name{Rule3}
\alias{Rule3}
\title{Rule3: Mean of survey data only}
\usage{
Rule3(data)
}
\arguments{
\item{data}{data.frame input called in OptRunCI (and subsequently used in DataGen)}
}
\value{
returns a mean value for the data.frame input 'Survey' column if data exists,
else returns NA
}
\description{
Rule3: Mean of survey data only
}
\details{
Internal GSO function. Rules functions appropriately combine observational
and expert data. These rules are used in the functions that generates scores
to use in the Growth Stage Optimisation (GSO).

Rule 3: This rule uses the mean of where the observational data exists,
otherwise uses nothing
}

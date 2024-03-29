% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSOFunc_gso.R
\name{gso}
\alias{gso}
\title{Growth Stage Optimisation (GSO)}
\usage{
gso(spp)
}
\arguments{
\item{spp}{the species lists from DataGen() function}
}
\value{
returns a data.frame with the result of the optimisation with a 'geom' column
}
\description{
Growth Stage Optimisation (GSO)
}
\details{
Internal GSO function.

Calcuates the GSO using a non-linear optimisation funciton. The input is the
species lists from DataGen() function
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fhProcess1.R
\name{fhProcess1}
\alias{fhProcess1}
\title{FIRE History Processing for FAME stage1}
\usage{
fhProcess1(
  inFH,
  inFHLayer = NULL,
  OtherAndUnknown = 2,
  validFIRETYPE = c("BURN", "BUSHFIRE", "UNKNOWN", "OTHER"),
  secondFH = NULL,
  secondFHLayer = NULL,
  baseFire = NULL,
  precision = 1
)
}
\arguments{
\item{inFH}{Input fire history polygon/mutipolygon data set with columns
SEASON and FIRETYPE  provided as shapefile, geopackage or ESRI geodatabase
file, or the name of an sf  object. if the geometries contain other than
polygon/mutipolygon an error will result.}

\item{inFHLayer}{Layer name if inFH has more than one layer
( for instance in a .gpkg) this allows selection of a particular layer,
otherwise fist layer is used (Default = NULL)}

\item{OtherAndUnknown}{integer Value to use for cases where fire type is:
"OTHER" or "UNKNOWN" = NA,
"BURN" = 1,
"BUSHFIRE" = 2.
NA = Fire excluded from analysis.
(usually set in settings file)}

\item{validFIRETYPE}{vector of valid names in the input FIRETYPE column in
the input fire history dataset(s), if the column contains NA or values not on
this list an error will occur}

\item{secondFH}{Second fire history to be combined with FH1 to make a fire
scenario same formats as for inFH}

\item{secondFHLayer}{Layer name if secondFH  has more than one layer ( for
instance in a .gpkg) this allows selection of a particular layer, otherwise
fist layer is used (Default = NULL)}

\item{baseFire}{Default NULL otherwise four digit integer SEASON  for fire
applied #' across the whole bounding box}

\item{precision}{default 10000, precision value fed through to GEOS see
\link[sf:st_precision]{sf::st_set_precision}, precision in metres = 1/precision 0 used precision
of exiting files}
}
\value{
inFH1 an sf geometry collection of processed fire history for input
into fhProcess2()
}
\description{
FIRE History Processing for FAME stage1
}
\details{
Takes a spatial polygon data file or r sf polygon dataset
checks that it contains the correct fields, and projection for FAME analysis
( using helper function fhCheck) including fields FIRETYPE and SEASON
Where polygons of different FIRETYPE or SEASON overlap the mainFHProcess()
function  called in this function constructs unique
non-overlapping polygon of their intersections ( and non intersecting areas
) and attributes each polygon with sequential fire SEASON (SEAS01, SEAS02
...) and corresponding FIRETYPE (TYPE01,TYPE02 ...)
}
\examples{
# randomFH<-generate_random_fire_history(20)
# plot(randomFH)
# outFH1<-fhProcess1(randomFH)
# plot(outFH1)
# #complete all fields in FH analysis
# outFH2<-fhProcess2(outFH1)
#
# plot(outFH2$OutDF,max.plot = 20)


}
\seealso{
\code{\link[=fhCheck]{fhCheck()}} for checking function run from within fhprocess1()
\link{mainFHProcess}  main processing function to derive Fire history analysis for
FAME
\link{prepFH} helper function that splits larger inFH files into gridded subunits
for parallel processing prior to passing to \link{mainFHProcess}
}

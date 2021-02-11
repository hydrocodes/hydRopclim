\name{indexcorrl}
\alias{indexcorrl}
\title{indexcorrl: Hydroclimatic indexes running correlation}
\usage{
indexcorrl(index.seas,variable.seas,rwin)
}
\description{
Running correlation between a climate index versus hydrological indexes, e.g., an ENSO index versus precipitation (Bourrel et al, 2015)
}
\arguments{
index.seas: A vector object, a yearly climate index

variable.seas: A matrix object, a yearly hydrological indexes

rwin: A numeric value, sliding window length in years
}
\examples{
indexcorrl(index.seas=p1, variable.seas=p2, rwin=11)
}

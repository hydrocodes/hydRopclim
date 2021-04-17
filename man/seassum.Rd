\name{seassum}
\alias{seassum}
\title{seassum: Calculation of seasonal sum matrix indexes for a season of n-months}
\usage{
seassum(p,start,win)
}
\description{
Calculation of seasonal sum matrix indexes for a season of n-months from m-hydroclimatic variables
}
\arguments{
p: A matrix object, monthly hydroclimate indexes

start: A numeric value, season starting month (e.g. September=9)

win: A numeric value, season length in months (e.g. 6 defines September-October-November-December-January-February)
}
\examples{
seassum(p=values3, start=9, win=6)
}

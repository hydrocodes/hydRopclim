\name{seasavg2}
\alias{seasavg2}
\title{seasavg2: Calculation of seasonal average matrix indexes for a season of n-months}
\usage{
seasavg2(p,start,win)
}
\description{
Calculation of seasonal average matrix indexes for a season of n-months from m-hydroclimatic variables
}
\arguments{
p: A matrix object, monthly hydroclimate indexes

start: A numeric value, season starting month (e.g. September=9)

win: A numeric value, season length in months (e.g. 3 defines September-October-November)
}
\examples{
seasavg2(p=values2, start=12, win=3)
}

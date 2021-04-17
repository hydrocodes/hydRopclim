\name{seasavg}
\alias{seasavg}
\title{seasavg: Calculation of seasonal average vector index for a season of n-months}
\usage{
seasavg(p,start,win)
}
\description{
Calculation of seasonal average vector index for a season of n-months
}
\arguments{
p: A vector object, a monthly climate index

start: A numeric value, season starting month (e.g. September=9)

win: A numeric value, season length in months (e.g. 3 defines September-October-November)
}
\examples{
seasavg(p=values1, start=9, win=3)
}

\name{tgridcorr}
\alias{tgridcorr}
\title{tgridcorr: Correcting temperature grid}
\usage{
tgridcorr(data,hBS,hNNR,hx,LR)
}
\description{
Correction of monthly grid-point mean temperature at target elevations, e.g., NCEP NCAR or others reanalysis products in function of an in-situ station (Rau et al, 2013).
}
\arguments{
data: A monthly dataframe object including Dates (Date), in-situ data (Station) and grid data at a level pressure (Grid.level)

hBS: A numeric value: station elevation in masl

hNNR: A numeric value: grid level elevation in masl

hx: A numeric value: Location elevation to correct in masl

LR: A 12 elements vector: mean monthly lapse rate in C/100m

}
\examples{
tgridcorr(data=database2, hBS=2479, hNNR=3012, hx=2000,
LR= c(-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6))
}
